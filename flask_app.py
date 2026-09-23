from __future__ import annotations

from flask import Flask, Response, render_template, request

app = Flask(__name__)

DEFAULTS = {
    "snp2genegrch38": "0",
    "chrcol": "NA",
    "poscol": "NA",
    "rsIDcol": "NA",
    "pcol": "NA",
    "eacol": "NA",
    "neacol": "NA",
    "orcol": "NA",
    "becol": "NA",
    "secol": "NA",
    "leadSNPsfile": "NA",
    "addleadSNPs": "1",
    "regionsfile": "NA",
    "N": "NA",
    "Ncol": "NA",
    "exMHC": "1",
    "MHCopt": "annot",
    "extMHC": "NA",
    "ensembl": "v102",
    "genetype": "protein_coding",
    "leadP": "5e-8",
    "gwasP": "0.05",
    "r2": "0.6",
    "r2_2": "0.1",
    "refpanel": "1KG/Phase3",
    "pop": "EUR",
    "MAF": "0",
    "refSNPs": "1",
    "mergeDist": "250",
    "magma": "0",
    "magma_window": "0",
    "magma_exp": "GTEx_v8/Brain_Amygdala.txt.gz:GTEx_v8/Brain_Anterior_cingulate_cortex_BA24.txt.gz",
    "posMap": "1",
    "posMapWindowSize": "NA",
    "posMapAnnot": "NA",
    "posMapCADDth": "0",
    "posMapRDBth": "NA",
    "posMapChr15": "NA",
    "posMapChr15Max": "NA",
    "posMapChr15Meth": "NA",
    "posMapAnnoDs": "NA",
    "posMapAnnoMeth": "NA",
}


def normalize_value(value: str | None, default: str = "NA") -> str:
    raw = (value or "").strip()
    return raw if raw else default


def build_params_config(form_data: dict | None = None) -> str:
    data = {**DEFAULTS}
    if form_data:
        for key, value in form_data.items():
            if key in data:
                data[key] = normalize_value(value, data[key])

    if data["snp2genegrch38"] not in {"0", "1"}:
        raise ValueError("snp2genegrch38 must be either 0 or 1")

    if data["N"] == "NA" and data["Ncol"] == "NA":
        data["N"] = "NA"
        data["Ncol"] = "NA"

    magma_exp = data["magma_exp"]
    if magma_exp:
        cleaned = ":".join(part.strip() for part in magma_exp.split(":") if part.strip())
        data["magma_exp"] = cleaned

    lines = [
        "[inputfiles]",
        f"snp2genegrch38 = {data['snp2genegrch38']}",
        f"chrcol = {data['chrcol']}",
        f"poscol = {data['poscol']}",
        f"rsIDcol = {data['rsIDcol']}",
        f"pcol = {data['pcol']}",
        f"eacol = {data['eacol']}",
        f"neacol = {data['neacol']}",
        f"orcol = {data['orcol']}",
        f"becol = {data['becol']}",
        f"secol = {data['secol']}",
        f"leadSNPsfile = {data['leadSNPsfile']}",
        f"addleadSNPs = {data['addleadSNPs']}",
        f"regionsfile = {data['regionsfile']}",
        "",
        "[params]",
        f"N = {data['N']}",
        f"Ncol = {data['Ncol']}",
        f"exMHC = {data['exMHC']}",
        f"MHCopt = {data['MHCopt']}",
        f"extMHC = {data['extMHC']}",
        f"ensembl = {data['ensembl']}",
        f"genetype = {data['genetype']}",
        f"leadP = {data['leadP']}",
        f"gwasP = {data['gwasP']}",
        f"r2 = {data['r2']}",
        f"r2_2 = {data['r2_2']}",
        f"refpanel = {data['refpanel']}",
        f"pop = {data['pop']}",
        f"MAF = {data['MAF']}",
        f"refSNPs = {data['refSNPs']}",
        f"mergeDist = {data['mergeDist']}",
        "",
        "[magma]",
        f"magma = {data['magma']}",
        f"magma_window = {data['magma_window']}",
        f"magma_exp = {data['magma_exp']}",
        "",
        "[posMap]",
        f"posMap = {data['posMap']}",
        f"posMapWindowSize = {data['posMapWindowSize']}",
        f"posMapAnnot = {data['posMapAnnot']}",
        f"posMapCADDth = {data['posMapCADDth']}",
        f"posMapRDBth = {data['posMapRDBth']}",
        f"posMapChr15 = {data['posMapChr15']}",
        f"posMapChr15Max = {data['posMapChr15Max']}",
        f"posMapChr15Meth = {data['posMapChr15Meth']}",
        f"posMapAnnoDs = {data['posMapAnnoDs']}",
        f"posMapAnnoMeth = {data['posMapAnnoMeth']}",
    ]
    return "\n".join(lines) + "\n"


@app.route("/", methods=["GET", "POST"])
def index():
    values = DEFAULTS.copy()
    error = None
    config_text = build_params_config(values)

    if request.method == "POST":
        values = {key: request.form.get(key, DEFAULTS.get(key, "NA")) for key in DEFAULTS}
        try:
            config_text = build_params_config(values)
        except ValueError as exc:
            error = str(exc)
            config_text = ""

    return render_template("index.html", values=values, config_text=config_text, error=error)


@app.route("/download", methods=["POST"])
def download():
    values = {key: request.form.get(key, DEFAULTS.get(key, "NA")) for key in DEFAULTS}
    try:
        config_text = build_params_config(values)
    except ValueError as exc:
        return Response(str(exc), status=400, mimetype="text/plain")

    return Response(
        config_text,
        mimetype="text/plain",
        headers={"Content-Disposition": "attachment; filename=params.config"},
    )


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
