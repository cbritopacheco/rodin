import os, subprocess, sys, json
EQS = {
 "eq1": r"\mathcal{L}P = -R, \qquad \mathcal{L}T = R - R_I, \qquad \mathcal{L}A = -R_I, \qquad \mathcal{L}G = -R_F",
 "eq2": r"R = k_a\,\frac{T^2}{K_A+T}\,\frac{P}{P_{\mathrm{ref}}}, \qquad R_I = k_{on} A T, \qquad R_F = k_F\,\frac{TG}{K_G+G}",
 "eq3": r"\mathcal{L}X = -k_{dX}X, \qquad -D\,\partial_n X = j_X E \ \ \text{(endocardio)}, \qquad \partial_t F = R_F",
 "eq4": r"R = \left[k_a\,\theta\,\frac{T^2}{K_A+T} + \kappa\,X\,\theta_X\right]\frac{P}{P_{\mathrm{ref}}}, \qquad \theta_X = \varepsilon + (1-\varepsilon)\frac{T}{K_A+T}",
 "eq5": r"\int_\Omega \frac{k_{TM}}{\delta}\,\chi\,(\tilde L - \tilde E)\,dV \;=\; k_{TM}\int_{\Gamma_w}(1-E)\,dA",
 "eq6": r"\Delta_{TM} = \min\!\left(A_n,\ T_n\frac{\Delta t\,k_{TM}^{\mathrm{loc}}}{1+\Delta t\,k_{TM}^{\mathrm{loc}}}\right), \qquad T \leftarrow T_n - \Delta_{TM}, \quad A \leftarrow A_n - \Delta_{TM}",
 "eq7": r"\frac{k_a}{k_{on}A_0} = \frac{0{,}035}{0{,}020} = 1{,}76",
}
TPL = r"""\documentclass[border=3pt,preview,12pt]{standalone}
\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{amsmath,amssymb}
\begin{document}$\displaystyle %s$\end{document}"""
out={}
for k,v in EQS.items():
    open(k+".tex","w").write(TPL % v)
    subprocess.run(["pdflatex","-interaction=nonstopmode",k+".tex"],capture_output=True)
    if not os.path.exists(k+".pdf"): print("FALLO",k); continue
    info=subprocess.run(["pdfinfo",k+".pdf"],capture_output=True,text=True).stdout
    w=[l for l in info.split("\n") if l.startswith("Page size")][0]
    wpt=float(w.split()[2])
    subprocess.run(["pdftoppm","-png","-r","300","-singlefile",k+".pdf",k],capture_output=True)
    out[k]=round(wpt/72*2.54,2)   # cm
json.dump(out,open("widths.json","w"))
print(out)
