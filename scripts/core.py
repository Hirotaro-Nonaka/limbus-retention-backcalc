"""Limbus retention back-calculation — reimplementation from Limbus_Retention_Reanalysis.md.
Shared core for tasks (a) pulse-shape sensitivity, (b) playtime-contamination bound, (c) 7-title generalization.
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy.optimize import nnls, minimize

DATA_DIR = Path(__file__).resolve().parent.parent / 'data' / 'steamdb'

# ---------- data prep ----------
def load_steamdb_daily(path):
    df = pd.read_csv(path, encoding='utf-8-sig')
    df['DateTime'] = pd.to_datetime(df['DateTime'])
    df['date'] = df['DateTime'].dt.date
    df['Average Players'] = pd.to_numeric(df['Average Players'], errors='coerce')
    d = df.groupby('date')['Average Players'].mean().dropna()
    d.index = pd.to_datetime(d.index)
    d = d[d > 0]
    # drop first day and (possibly incomplete) last day
    d = d.iloc[1:-1]
    return d

def weekday_adjust(s):
    ma = s.rolling(7, center=True, min_periods=4).mean()
    ratio = s / ma
    coef = ratio.groupby(ratio.index.dayofweek).median()
    adj = s / s.index.dayofweek.map(coef).values
    return adj, coef

# ---------- kernels ----------
def R_kernel(T, lam, k, p, c, tau):
    d = np.arange(T, dtype=float)
    r = (1-p)*np.exp(-(d/lam)**k) + p*c*np.exp(-d/tau)
    r[0] = 1.0
    return r

def A_kernel(T, lam, k, p, tau):
    d = np.arange(T, dtype=float)
    x = np.maximum(0.0, d-29)
    return (1-p)*np.exp(-(x/lam)**k) + p*np.exp(-x/tau)

# ---------- pulse shapes ----------
def pulse_shape(name, W):
    j = np.arange(W, dtype=float)
    if name == 'exp3':   s = np.exp(-j/3)
    elif name == 'exp2': s = np.exp(-j/2)
    elif name == 'exp5': s = np.exp(-j/5)
    elif name == 'rect': s = np.ones(W)
    elif name == 'tri':  s = 1 - j/W
    elif name == 'gamma':s = (j+1)*np.exp(-j/2)   # humped, peak day1
    elif name == 'delta':s = np.zeros(W); s[0] = 1
    else: raise ValueError(name)
    return s/s.sum()   # normalize so amplitude = total activations

def build_X(T, events, R, shape='exp3', W_default=7, W_long=14):
    """events: list of (day_index, is_long). Returns T x (n_ev+1) matrix; col0 = Base(const inflow)."""
    cols = [np.cumsum(R)]  # base: unit inflow/day convolved with R
    for t0, is_long in events:
        W = W_long if is_long else W_default
        if shape == 'delta': W = 1
        s = pulse_shape(shape, W)
        col = np.zeros(T)
        for j, sj in enumerate(s):
            if t0+j >= T: break
            col[t0+j:] += sj * R[:T-t0-j]
        cols.append(col)
    return np.column_stack(cols)

def build_X_kernel(T, events, K, shape='exp3', W_default=7, W_long=14):
    """same but with arbitrary kernel K (e.g. MAU kernel)."""
    cols = [np.cumsum(K)]
    for t0, is_long in events:
        W = W_long if is_long else W_default
        if shape == 'delta': W = 1
        s = pulse_shape(shape, W)
        col = np.zeros(T)
        for j, sj in enumerate(s):
            if t0+j >= T: break
            col[t0+j:] += sj * K[:T-t0-j]
        cols.append(col)
    return np.column_stack(cols)

# ---------- objective ----------
class Fitter:
    def __init__(self, y, events, shape='exp3', anchors=None, beta=0.0,
                 h_prior=(2.2, 0.2), use_h=True, W=(7,14)):
        """y: adjusted daily avg CCU (np array). events: [(idx,is_long)].
        anchors: dict(jan_mask=bool array, dau=..., mau=..., ratio=(0.455,0.02)) or None.
        beta: playtime-contamination coefficient (task b): CCU = DAU*h/24*(1+beta*g(t)),
              g = unit-peak pulse profile summed over events."""
        self.y = y; self.T = len(y); self.events = events; self.shape = shape
        self.anchors = anchors; self.beta = beta; self.h_prior = h_prior; self.use_h = use_h; self.W = W
        g = np.zeros(self.T)
        for t0, is_long in events:
            W = 14 if is_long else 7
            s = pulse_shape('exp3', W); s = s/s[0]  # unit peak
            for j, sj in enumerate(s):
                if t0+j < self.T: g[t0+j] = max(g[t0+j], sj)
        self.g = g

    def unpack(self, th):
        lam = np.exp(th[0]); k = th[1]; p = th[2]; c = th[3]; tau = np.exp(th[4])
        h = th[5] if self.use_h else 2.2
        return lam, k, p, c, tau, h

    def solve(self, th):
        lam, k, p, c, tau, h = self.unpack(th)
        R = R_kernel(self.T, lam, k, p, c, tau)
        X = build_X(self.T, self.events, R, self.shape, self.W[0], self.W[1])
        mult = (h/24.0) * (1.0 + self.beta*self.g)          # CCU = DAU * h/24 * (1+beta g)
        rows = [X * (mult/self.y)[:, None]]
        rhs  = [np.ones(self.T)]
        extra = 0.0
        if self.anchors:
            A = A_kernel(self.T, lam, k, p, tau)
            Xm = build_X_kernel(self.T, self.events, A, self.shape, self.W[0], self.W[1])
            m = self.anchors['jan_mask']
            dau_row = X[m].mean(axis=0); mau_row = Xm[m].mean(axis=0)
            rows.append((dau_row/(0.05*self.anchors['dau']))[None, :]); rhs.append(np.array([1/0.05]))
            rows.append((mau_row/(0.03*self.anchors['mau']))[None, :]); rhs.append(np.array([1/0.03]))
        Xw = np.vstack(rows); yw = np.concatenate(rhs)
        a, rnorm = nnls(Xw, yw, maxiter=10*Xw.shape[1])
        loss = rnorm**2
        if self.anchors:
            A = A_kernel(self.T, lam, k, p, tau)
            Xm = build_X_kernel(self.T, self.events, A, self.shape, self.W[0], self.W[1])
            m = self.anchors['jan_mask']
            dau = X[m].mean(axis=0) @ a; mau = Xm[m].mean(axis=0) @ a
            r0, rs = self.anchors.get('ratio', (0.455, 0.02))
            loss += ((dau/mau - r0)/rs)**2
        if self.use_h:
            loss += ((h - self.h_prior[0])/self.h_prior[1])**2
        return loss, a

    def obj(self, th):
        try:
            return self.solve(th)[0]
        except Exception:
            return 1e12

    BOUNDS = [(np.log(0.1), np.log(30)), (0.08, 1.5), (1e-3, 0.4), (0.2, 1.0),
              (np.log(200), np.log(6000)), (1.0, 4.0)]

    def obj_b(self, th):
        lo = np.array([b[0] for b in self.BOUNDS]); hi = np.array([b[1] for b in self.BOUNDS])
        thc = np.clip(th, lo, hi)
        pen = 1e3*np.sum((th-thc)**2)
        try:
            return self.solve(thc)[0] + pen
        except Exception:
            return 1e12

    def fit(self, th0=None, maxiter=3000):
        if th0 is None:
            th0 = np.array([np.log(1.5), 0.30, 0.023, 0.90, np.log(6000), 2.2])
        res = minimize(self.obj_b, th0, method='Powell', bounds=self.BOUNDS,
                       options=dict(maxiter=maxiter, xtol=1e-5, ftol=1e-7))
        th = np.clip(res.x, [b[0] for b in self.BOUNDS], [b[1] for b in self.BOUNDS])
        res2 = minimize(self.obj_b, th, method='Nelder-Mead',
                        options=dict(maxiter=600, xatol=1e-5, fatol=1e-8))
        th2 = np.clip(res2.x, [b[0] for b in self.BOUNDS], [b[1] for b in self.BOUNDS])
        if self.obj_b(th2) < self.obj_b(th):
            th = th2
        return th, res

    def report(self, th):
        lam, k, p, c, tau, h = self.unpack(th)
        loss, a = self.solve(th)
        R = R_kernel(max(self.T, 400), lam, k, p, c, tau)
        X = build_X(self.T, self.events, R_kernel(self.T, lam, k, p, c, tau), self.shape, self.W[0], self.W[1])
        mult = (h/24.0)*(1.0+self.beta*self.g)
        pred = (X @ a) * mult
        rmse = np.sqrt(np.mean(((pred-self.y)/self.y)**2))
        out = dict(D1=R[1], D7=R[7], D30=R[30], D180=R[180], D365=R[365] if len(R) > 365 else np.nan,
                   lam=lam, k=k, p=p, c=c, tau=tau, h=h, rmse=rmse, loss=loss,
                   base=a[0], amps=a[1:], pred=pred)
        if self.anchors:
            A = A_kernel(self.T, lam, k, p, tau)
            Xm = build_X_kernel(self.T, self.events, A, self.shape, self.W[0], self.W[1])
            m = self.anchors['jan_mask']
            out['jan_dau'] = float(X[m].mean(axis=0) @ a)
            out['jan_mau'] = float(Xm[m].mean(axis=0) @ a)
        return out

# ---------- Limbus event calendar (from Reanalysis md + inflow.csv spike dates) ----------
LIMBUS_EVENTS = {
 'Launch': ('2023-02-27', True),
 'S2/CantoIV': ('2023-06-01', False), 'S3/CantoV': ('2023-11-16', False),
 'S4/CantoVI': ('2024-03-28', False), 'S5/CantoVII': ('2024-10-10', False),
 'S6/CantoVIII': ('2025-05-15', False), 'S7/CantoIX': ('2025-12-31', False),
 'S2/Story': ('2023-06-15', False), 'S3/Story': ('2023-11-30', False),
 'S4/Story': ('2024-04-11', False), 'S5/Story': ('2024-10-24', False),
 'S6/Story': ('2025-05-29', False), 'S7/Story': ('2026-01-14', False),
 'WN1': ('2023-10-26', False), 'WN2': ('2024-01-11', False), 'WN3': ('2024-05-02', False),
 'WN4': ('2024-09-05', False), 'WN5': ('2025-01-09', False), 'WN6': ('2025-07-17', False),
 'WN7': ('2025-11-06', False), 'WN8': ('2026-03-05', False),
 'Collab': ('2025-09-25', True),
}
SPIKES = ['2023-03-23','2023-04-06','2023-04-20','2023-06-09','2023-06-29','2023-07-13','2023-07-27',
'2023-08-24','2023-09-07','2023-09-15','2023-09-23','2023-11-24','2023-12-14','2023-12-22','2023-12-30',
'2024-02-01','2024-02-22','2024-04-05','2024-05-16','2024-06-13','2024-06-27','2024-07-25','2024-08-08',
'2024-08-22','2024-10-03','2024-10-19','2024-11-14','2024-11-23','2024-12-12','2024-12-26','2025-01-23',
'2025-02-06','2025-02-20','2025-02-28','2025-03-08','2025-03-20','2025-05-23','2025-07-31','2025-08-14',
'2025-08-28','2026-01-08','2026-02-19','2026-03-13','2026-04-16','2026-04-30','2026-05-14','2026-05-28',
'2026-06-11','2026-06-25']

def limbus_setup(shape='exp3', beta=0.0, W=(7,14)):
    s = load_steamdb_daily(DATA_DIR / 'Limbus_steamdb_chart_1973530.csv')
    adj, _ = weekday_adjust(s)
    idx = {d: i for i, d in enumerate(adj.index)}
    events = []
    for name, (d, lng) in LIMBUS_EVENTS.items():
        events.append((idx[pd.Timestamp(d)], lng))
    for d in SPIKES:
        events.append((idx[pd.Timestamp(d)], False))
    jan = (adj.index >= '2026-01-01') & (adj.index <= '2026-01-31')
    anchors = dict(jan_mask=np.asarray(jan), dau=360639.0, mau=792250.0, ratio=(0.455, 0.02))
    f = Fitter(adj.values.astype(float), events, shape=shape, anchors=anchors, beta=beta, W=W)
    return f, adj
