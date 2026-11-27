# Column specifications and header rows of the tables, copied verbatim from
# the LaTeX sources (generated once; edit here to change a table's header).
paper_tables <- list(
  tab_E7 = list(spec = r"---[rrrcccccccc]---",
    header = r"---[$S$ & $\varepsilon$ & $\neff$ & Kolmogorov & $0.04/\sqrt{\neff}$ & $3.1/\sqrt{\neff}$ & $\log_{10}\PP(\Zhat\in\{0,\infty\})$ & $\log_{10}e^{-4N\Hov}$ & coverage & $\PP(|\Delta|\ge0.25)$ & bound\\]---"),
  tab_EB = list(spec = r"---[lrrrcccc]---",
    header = r"---[Family & $S$ & med $\neffhat$ & med $\Ihat$ & coverage (delta) & coverage, \cref{prop:ebci} & bounded & width ratio\\]---"),
  tab_NP = list(spec = r"---[rrccccccccccc]---",
    header = r"---[$d$ & $S$ & $\beta$ (pred.) & med $1-\Ihat$ & $\neff$ (pred.) & med $\neffhat$ & RMSE & med MCSE & coverage & flag (\cref{alg:diag}) & half-width, \cref{prop:ebci} & coverage, \cref{prop:ebci} & med $\widehat\Lambda_3$\\]---"),
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
  tab_prop = list(spec = r"---[rlccccrr]---",
    header = r"---[$d$ & proposal & sd of $\log\Zhat$ & bias/MCSE & RMSE & coverage & med $\neffhat$ & $t_\Tp$ rejections\\]---"),
  tab_E5 = list(spec = r"---[lrccccc]---",
    header = r"---[Scenario & med $\neffhat$ & bias/MCSE & Var/MCSE$^2$ & coverage & mean $\Vjack$ & variance check (\%)\\]---"),
)