# Column specifications and header rows of the tables, copied verbatim from
# the LaTeX sources (generated once; edit here to change a table's header).
paper_tables <- list(
  tab_ES_top = list(spec = r"---[llrcccrcrr]---",
    header = r"---[Parameterization & draws & $S$ & coverage & RMSE & med MCSE & med $\neffhat$ & med $\Ihat$ & overlap flag & failure flags\\]---"),
  tab_ES_bottom = list(spec = r"---[rlrcrrc]---",
    header = r"---[$S$ & estimator of $\tau$ & med $\widehat\tau$ & med MCSE & med $\neffPhat$ & flag (\%) & coverage\\]---"),
  tab_E1_top = list(spec = r"---[lccccccc]---",
    header = r"---[Problem & iter.\ ($\pm20$) & iter.\ ($\pm200$) & final error & $\sup_x|\varphi'(x)|$ & $T'$ (formula) & $T'$ (fin.\ diff.) & $1-\Hhat/\Hzero$\\]---"),
  tab_E1_bottom = list(spec = r"---[rrrrr]---",
    header = r"---[$\Sone=\Stwo$ & $\neff$ & mean $T'$ & sd $T'$ & $\sqrt{\neff}\,$sd $T'$\\]---"),
  tab_E7 = list(spec = r"---[rrrcccccccc]---",
    header = r"---[$S$ & $\varepsilon$ & $\neff$ & Kolmogorov & $0.04/\sqrt{\neff}$ & $3.1/\sqrt{\neff}$ & $\log_{10}\PP(\Zhat\in\{0,\infty\})$ & $\log_{10}e^{-4N\Hov}$ & coverage & $\PP(|\Delta|\ge0.25)$ & bound\\]---"),
  tab_E6 = list(spec = r"---[llrrrccccc]---",
    header = r"---[Family & parameters & $\Sone=\Stwo$ & $\Hov$ & $\neff$ & coverage & KS & $\mathrm{sd}(\Delta)/\sigN$ & med $|\log\frac{\sighat}{\sigN}|$ & med $\neffhat/\neff$\\]---"),
  tab_EB = list(spec = r"---[lrrrcccc]---",
    header = r"---[Family & $S$ & med $\neffhat$ & med $\Ihat$ & coverage (delta) & coverage, \cref{prop:ebci} & bounded & width ratio\\]---"),
  tab_NP = list(spec = r"---[rrccccccccccc]---",
    header = r"---[$d$ & $S$ & $\beta$ (pred.) & med $1-\Ihat$ & $\neff$ (pred.) & med $\neffhat$ & RMSE & med MCSE & coverage & flag (\cref{alg:diag}) & half-width, \cref{prop:ebci} & coverage, \cref{prop:ebci} & med $\widehat\Lambda_3$\\]---"),
  tab_ADV = list(spec = r"---[lrrrccc]---",
    header = r"---[Example & $S$ & $\neff$ & med $\neffhat$ & med $\Ihat$ & coverage & KS\\]---"),
  tab_E2 = list(spec = r"---[lrrrcccccccccc]---",
    header = r"---[Setting & $S$ & $G$ & med $\neffhat$ & Var/MCSE$^2$ & mean $\Vjack$ & $\Vjack$ & $\Imax$ (Bonf.) & $\Imax$ (naive) & $\Imax$ (exact) & $\Gr$ & $t_\Tp$ & coverage\\]---"),
  tab_FA = list(spec = r"---[lccccccc]---",
    header = r"---[Setting & $t_\Tp$ & $\Vjack$ & $\Imax$ & $\Gr$ & jackknife union & union of four & union, $\alpha/4$ each\\]---"),
  tab_E3 = list(spec = r"---[c|ccccccc|cccc]---",
    header = r"---[& \multicolumn{7}{c|}{spectral MCSE} & \multicolumn{4}{c}{naive MCSE}\\
$\rho$ & Var/MCSE$^2$ & coverage & mean $\Vjack$ & $\Vjack$ & $\Imax$ & $\Gr$ & $t_\Tp$ & Var/MCSE$^2$ & coverage & mean $\Vjack$ & $\Vjack$\\]---"),
  tab_E8 = list(spec = r"---[rrrrrrcccc]---",
    header = r"---[$\gamma$ & $S$ & $\neff$ & $\neff^P$ & $\log_{10}\PP(\Zhat\in\{0,\infty\})$ & lower bound & Kolmogorov & $\gamma\sqrt{v_P}\times$Kolm. & coverage ($\sigma_{N,P}$) & coverage ($\sigN$)\\]---"),
  tab_E9_top = list(spec = r"---[llrccccccccc]---",
    header = r"---[& & & \multicolumn{6}{c}{coverage, spectral MCSE, at increasing $\neffP$} & & \\
\cmidrule(lr){4-9}
Family & $\rho$ & $\tau_W$ & $1$ & $2$ & $4$ & $8$ & $20$ & $80$ & oracle & ignoring autocorr.\\]---"),
  tab_E9_bottom = list(spec = r"---[llrccccccccc]---",
    header = r"---[& & & \multicolumn{6}{c}{median $\widehat n_P/v_P$ at increasing $\neffP$} & \multicolumn{3}{c}{median $\neffPhat$}\\
\cmidrule(lr){4-9}\cmidrule(lr){10-12}
Family & $\rho$ & $\tau_W$ & $1$ & $2$ & $4$ & $8$ & $20$ & $80$ & $1$ & $2$ & $4$\\]---"),
  tab_LRV = list(spec = r"---[lrrrrccccrr]---",
    header = r"---[Family & $\rho$ & $\tau_W$ & $\neffP$ & $S$ & oracle & Geyer & AR & \texttt{error\_measures} & med $\neffPhat$ (Geyer) & med $\neffPhat$ (AR)\\]---"),
  tab_E4 = list(spec = r"---[rrrcccccccccc]---",
    header = r"---[& & & \multicolumn{3}{c}{sd of $\log\Zhat$} & \multicolumn{2}{c}{bias / MCSE} & \multicolumn{2}{c}{$R_g$} & \multicolumn{3}{c}{mean $\Vjack$}\\
\cmidrule(lr){4-6}\cmidrule(lr){7-8}\cmidrule(lr){9-10}\cmidrule(lr){11-13}
$d$ & $S$ & reps & fixed & same & split & same & split & draws fixed & pathwise & no refit & refit, draws fixed & refit, CRN\\]---"),
  tab_prop = list(spec = r"---[rlccccrr]---",
    header = r"---[$d$ & proposal & sd of $\log\Zhat$ & bias/MCSE & RMSE & coverage & med $\neffhat$ & $t_\Tp$ rejections\\]---"),
  tab_E4c = list(spec = r"---[rrrcrcc]---",
    header = r"---[$d$ & $k$ & $S$ & bias/MCSE & med $\neffhat$ & med $\Ihat$ & $\widehat C$\\]---"),
  tab_E5 = list(spec = r"---[lrccccccccc]---",
    header = r"---[Scenario & med $\neffhat$ & bias/MCSE & Var/MCSE$^2$ & coverage & mean $\Vjack$ & $\Vjack$ & $\Imax$ & $\Gr$ & $t_\Tp$\\]---"),
  tab_R1 = list(spec = r"---[lrrrrcrrrr]---",
    header = r"---[Posterior and proposal & $\Sone$ & $\Stwo$ & $\neff$ & $N\Hov$ & max ratio & $\PP(\neffhat\ge10)$ & bound & $q_{0.99}(\neffhat)$ & bound\\]---"),
  tab_R2 = list(spec = r"---[lrrrcrr]---",
    header = r"---[Posterior and proposal & $\Sone$ & $\Stwo$ & $\neff$ & sd ratio & rejection & sd$(t_\Tp)$\\]---")
)
