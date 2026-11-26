import numpy as np
import pandas as pd
from itertools import combinations
from math import comb
from scipy.stats import kendalltau
from pathlib import Path
OUT=Path(__file__).resolve().parent
P=[1,2,3,5,10]; B=100; NS=[100,500]; SEED=20261126
SC=['independent','moderate','strong','duplicate','correlated_duplicate','nonlinear']
def norm(X): return (X-X.min(0))/np.ptp(X,axis=0)
def centered(X):
    V=X-X.mean(0); return V/np.linalg.norm(V,axis=0)
def residual(U,j,S):
    if not len(S): return 1.
    a,d,_=np.linalg.svd(U[:,S],full_matrices=False)
    Q=a[:,d>1e-10*d[0]]; r=U[:,j]-Q@(Q.T@U[:,j]); return np.clip(r@r,0,1)
def weights(Z):
    m=Z.shape[1]; U=centered(Z); q=np.ones(m)/m
    for k in range(1,m):
        a=1/(m*comb(m-1,k))
        for S in combinations(range(m),k):
            Q,d,_=np.linalg.svd(U[:,S],full_matrices=False); Q=Q[:,d>1e-10*d[0]]
            J=[j for j in range(m) if j not in S]; R=U[:,J]-Q@(Q.T@U[:,J]); q[J]+=a*np.clip(np.sum(R**2,axis=0),0,1)
    return q,q/q.sum()
def scores(Z,w): return np.column_stack([(Z**p@w)**(1/p) for p in P])
def dp2(X,classic=True,tol=1e-8,maxit=100):
    D=(X-X.min(0))/X.std(0,ddof=1) if classic else X
    U=centered(D); ci=D.sum(1); seen=set(); status='maxit'; m=D.shape[1]; q=np.ones(m)
    for k in range(maxit):
        a=ci-ci.mean(); rho=np.abs(U.T@a/np.linalg.norm(a)); order=tuple(np.argsort(-rho,kind='stable'))
        if order in seen: status='cycle'; break
        seen.add(order); q=np.ones(m)
        for t,j in enumerate(order): q[j]=residual(U,j,order[:t])
        new=D@q
        if abs(new.mean()-ci.mean())<=tol: ci=new; status='mean_tolerance'; break
        ci=new
        a=ci-ci.mean(); ordernew=tuple(np.argsort(-np.abs(U.T@a/np.linalg.norm(a)),kind='stable'))
        if ordernew==order: status='fixed_order'; break
    return ci,q/q.sum(),status
def kt(a,b): return kendalltau(a,b).statistic
rng=np.random.default_rng(SEED); weightrows=[]; metricrows=[]; duprows=[]; statusrows=[]
for n in NS:
    for b in range(B):
        a,u,c,d=rng.random((4,n)); indep=np.column_stack([a,u,c,d]); moderate=np.column_stack([a,.5*a+.5*u,c,d]); strong=np.column_stack([a,.9*a+.1*u,c,d]); datasets=[indep,moderate,strong,np.column_stack([indep,a]),np.column_stack([strong,a]),np.column_stack([a,(2*a-1)**2,c,d])]
        keep=np.sort(rng.choice(n,size=int(.9*n),replace=False)); cache={}
        for s,X in zip(SC,datasets):
            Z=norm(X); q,w=weights(Z); C=scores(Z,w); dc,wc,st=dp2(X); dm,wm,sm=dp2(Z,False); cache[s]=(C,dc,dm,w,wc,wm)
            statusrows.append([n,b,s,st,sm]); Zr=norm(X[keep]); _,wr=weights(Zr); Cr=scores(Zr,wr); dcr,_,strr=dp2(X[keep]); dmr,_,smr=dp2(Zr,False)
            for j in range(len(w)): weightrows.append([n,b,s,j+1,q[j],w[j],wc[j],wm[j]])
            for j,p in enumerate(P):
                metricrows.append([n,b,s,p,kt(C[:,0],C[:,j]),kt(C[:,j],dc),kt(C[:,j],dm),kt(C[keep,j],Cr[:,j]),kt(dc[keep],dcr),kt(dm[keep],dmr),np.mean(C[:,j]),np.mean(C[:,j]-C[:,0]),np.mean(np.abs(w-wr)),st,sm])
            if s in ('duplicate','correlated_duplicate'):
                base='independent' if s=='duplicate' else 'strong'; cb,db,mb,wb,wcb,wmb=cache[base]
                for j,p in enumerate(P): duprows.append([n,b,s,p,kt(cb[:,j],C[:,j]),np.max(np.abs(cb[:,j]-C[:,j])),w[0]+w[-1]-wb[0],kt(db,dc),kt(mb,dm),wc[0]+wc[-1]-wcb[0],wm[0]+wm[-1]-wmb[0]])
    print('completed',n,flush=True)
W=pd.DataFrame(weightrows,columns='n rep scenario indicator q weight dp2_q_norm dp2_matched_weight'.split()); M=pd.DataFrame(metricrows,columns='n rep scenario p tau_p1 tau_dp2 tau_dp2_matched tau_delete tau_delete_dp2 tau_delete_matched mean_score mean_gap mean_abs_weight_change status_dp2 status_matched'.split()); D=pd.DataFrame(duprows,columns='n rep scenario p tau_dup max_abs_score_change group_weight_change tau_dup_dp2 tau_dup_matched group_change_dp2_q group_change_matched'.split()); T=pd.DataFrame(statusrows,columns='n rep scenario status_dp2 status_matched'.split())
for df,name,keys in [(W,'weights',['n','scenario','indicator']),(M,'metrics',['n','scenario','p']),(D,'duplication',['n','scenario','p'])]:
    df.to_csv(OUT/(name+'_raw.csv'),index=False); nums=[c for c in df.select_dtypes('number').columns if c not in keys+['rep']]; g=df.groupby(keys,sort=False)[nums]; summary=g.mean().join(g.std().add_suffix('_sd')).join(g.quantile(.05).add_suffix('_q05')).join(g.quantile(.95).add_suffix('_q95')); summary.to_csv(OUT/(name+'_summary.csv'))
T.to_csv(OUT/'dp2_status.csv',index=False)
