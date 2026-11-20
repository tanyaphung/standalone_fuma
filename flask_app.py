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
    "magma_exp": 'GTEx/v8/gtex_v8_ts_avg_log2TPM:GTEx/v8/gtex_v8_ts_general_avg_log2TPM',
    "posMap": "1",
    "posMapWindowSize": "NA",
    "posMapAnnot": "NA",
    "posMapCADDcheck": "0",
    "posMapCADDth": "0",
    "posMapRDBcheck": "0",
    "posMapRDBth": "NA",
    "posMapChr15check": "0",
    "posMapChr15": "NA",
    "posMapChr15Max": "NA",
    "posMapChr15Meth": "NA",
    "posMapAnnoDs": "NA",
    "posMapAnnoMeth": "NA",
    "eqtlMap": "0",
    "eqtlMaptss": "NA",
    "eqtlMapSig": "1",
    "eqtlMapP": "1",
    "eqtlMapCADDth": "0",
    "eqtlMapRDBth": "NA",
    "eqtlMapChr15": "NA",
    "eqtlMapChr15Max": "NA",
    "eqtlMapChr15Meth": "NA",
    "eqtlMapAnnoDs": "NA",
    "eqtlMapAnnoMeth": "NA",
    "xqtlsMap": "0",
    "xqtlsMapdss": "NA",
    "xqtlP": "1e-3",
    "ciMap": "0",
    "ciMapBuiltin": "NA",
    "ciMapFileN": "0",
    "ciMapFiles": "NA",
    "ciMapFDR": "NA",
    "ciMapPromWindow": "NA",
    "ciMapRoadmap": "NA",
    "ciMapEnhFilt": "0",
    "ciMapPromFilt": "0",
    "ciMapCADDth": "0",
    "ciMapRDBth": "NA",
    "ciMapChr15": "NA",
    "ciMapChr15Max": "NA",
    "ciMapChr15Meth": "NA",
    "ciMapAnnoDs": "NA",
    "ciMapAnnoMeth": "NA",
}


def normalize_value(value: str | None, default: str = "NA") -> str:
    raw = (value or "").strip()
    return raw if raw else default


def normalize_multiselect(value, default: str = "NA") -> str:
    if value is None:
        return default
    if isinstance(value, (list, tuple, set)):
        parts = [str(v).strip() for v in value if str(v).strip()]
        return ":".join(parts) if parts else default
    raw = str(value).strip()
    return raw if raw else default


def build_params_config(form_data: dict | None = None) -> str:
    data = {**DEFAULTS}
    if form_data:
        for key, value in form_data.items():
            if key in data:
                if key in {"posMapAnnot", "xqtlsMapdss"}:
                    data[key] = normalize_multiselect(value, data[key])
                else:
                    data[key] = normalize_value(value, data[key])

    if data["posMapCADDcheck"] == "0":
        data["posMapCADDth"] = "0"
    if data["posMapRDBcheck"] == "0":
        data["posMapRDBth"] = "NA"
    elif data["posMapRDBth"] == "NA":
        data["posMapRDBth"] = "7"
    if data["posMapChr15check"] == "0":
        data["posMapChr15Max"] = "NA"
    elif data["posMapChr15Max"] == "NA":
        data["posMapChr15Max"] = "7"
    if data["posMapChr15check"] == "0":
        data["ciMapChr15Meth"] = "NA"
    elif data["ciMapChr15Meth"] == "NA":
        data["ciMapChr15Meth"] = "any"

    if data["snp2genegrch38"] not in {"0", "1"}:
        raise ValueError("snp2genegrch38 must be either 0 or 1")

    n_value = normalize_value(data.get("N"))
    ncol_value = normalize_value(data.get("Ncol"))
    if n_value != "NA":
        data["N"] = n_value
        data["Ncol"] = "NA"
    else:
        data["N"] = "NA"
        data["Ncol"] = ncol_value

    pos_map_window = normalize_value(data.get("posMapWindowSize"))
    pos_map_annot = normalize_value(data.get("posMapAnnot"))
    if pos_map_window != "NA":
        data["posMapWindowSize"] = pos_map_window
        data["posMapAnnot"] = "NA"
    else:
        data["posMapWindowSize"] = "NA"
        data["posMapAnnot"] = pos_map_annot

    magma_exp = data["magma_exp"]
    if magma_exp:
        cleaned = ":".join(part.strip() for part in magma_exp.split(":") if part.strip())
        data["magma_exp"] = cleaned

    lines = [
        "[inputfiles]",
        f"snp2genegrch38={data['snp2genegrch38']}",
        f"chrcol={data['chrcol']}",
        f"poscol={data['poscol']}",
        f"rsIDcol={data['rsIDcol']}",
        f"pcol={data['pcol']}",
        f"eacol={data['eacol']}",
        f"neacol={data['neacol']}",
        f"orcol={data['orcol']}",
        f"becol={data['becol']}",
        f"secol={data['secol']}",
        f"leadSNPsfile={data['leadSNPsfile']}",
        f"addleadSNPs={data['addleadSNPs']}",
        f"regionsfile={data['regionsfile']}",
        "",
        "[params]",
        f"N={data['N']}",
        f"Ncol={data['Ncol']}",
        f"exMHC={data['exMHC']}",
        f"MHCopt={data['MHCopt']}",
        f"extMHC={data['extMHC']}",
        f"ensembl={data['ensembl']}",
        f"genetype={data['genetype']}",
        f"leadP={data['leadP']}",
        f"gwasP={data['gwasP']}",
        f"r2={data['r2']}",
        f"r2_2={data['r2_2']}",
        f"refpanel={data['refpanel']}",
        f"pop={data['pop']}",
        f"MAF={data['MAF']}",
        f"refSNPs={data['refSNPs']}",
        f"mergeDist={data['mergeDist']}",
        "",
        "[magma]",
        f"magma={data['magma']}",
        f"magma_window={data['magma_window']}",
        f"magma_exp={data['magma_exp']}",
        "",
        "[posMap]",
        f"posMap={data['posMap']}",
        f"posMapWindowSize={data['posMapWindowSize']}",
        f"posMapAnnot={data['posMapAnnot']}",
        f"posMapCADDth={data['posMapCADDth']}",
        f"posMapRDBth={data['posMapRDBth']}",
        f"posMapChr15={data['posMapChr15']}",
        f"posMapChr15Max={data['posMapChr15Max']}",
        f"posMapChr15Meth={data['posMapChr15Meth']}",
        f"posMapAnnoDs={data['posMapAnnoDs']}",
        f"posMapAnnoMeth={data['posMapAnnoMeth']}",
        "",
        "[eqtlMap]",
        f"eqtlMap={data['eqtlMap']}",
        f"eqtlMaptss={data['eqtlMaptss']}",
        f"eqtlMapSig={data['eqtlMapSig']}",
        f"eqtlMapP={data['eqtlMapP']}",
        f"eqtlMapCADDth={data['eqtlMapCADDth']}",
        f"eqtlMapRDBth={data['eqtlMapRDBth']}",
        f"eqtlMapChr15={data['eqtlMapChr15']}",
        f"eqtlMapChr15Max={data['eqtlMapChr15Max']}",
        f"eqtlMapChr15Meth={data['eqtlMapChr15Meth']}",
        f"eqtlMapAnnoDs={data['eqtlMapAnnoDs']}",
        f"eqtlMapAnnoMeth={data['eqtlMapAnnoMeth']}",
        "",
        "[xqtlsMap]",
        f"xqtlsMap={data['xqtlsMap']}",
        f"xqtlsMapdss={data['xqtlsMapdss']}",
        f"xqtlP={data['xqtlP']}",
        "",
        "[ciMap]",
        f"ciMap={data['ciMap']}",
        f"ciMapBuiltin={data['ciMapBuiltin']}",
        f"ciMapFileN={data['ciMapFileN']}",
        f"ciMapFiles={data['ciMapFiles']}",
        f"ciMapFDR={data['ciMapFDR']}",
        f"ciMapPromWindow={data['ciMapPromWindow']}",
        f"ciMapRoadmap={data['ciMapRoadmap']}",
        f"ciMapEnhFilt={data['ciMapEnhFilt']}",
        f"ciMapPromFilt={data['ciMapPromFilt']}",
        f"ciMapCADDth={data['ciMapCADDth']}",
        f"ciMapRDBth={data['ciMapRDBth']}",
        f"ciMapChr15={data['ciMapChr15']}",
        f"ciMapChr15Max={data['ciMapChr15Max']}",
        f"ciMapChr15Meth={data['ciMapChr15Meth']}",
        f"ciMapAnnoDs={data['ciMapAnnoDs']}",
        f"ciMapAnnoMeth={data['ciMapAnnoMeth']}",
    ]
    return "\n".join(lines) + "\n"


@app.route("/", methods=["GET", "POST"])
def index():
    values = DEFAULTS.copy()
    error = None
    config_text = build_params_config(values)

    if request.method == "POST":
        values = DEFAULTS.copy()
        for key in DEFAULTS:
            if key == "addleadSNPs":
                values[key] = "1" if request.form.get("addleadSNPs") == "1" else "0"
            elif key in {"posMapAnnot", "xqtlsMapdss"}:
                values[key] = ":".join(request.form.getlist(key))
            else:
                values[key] = request.form.get(key, DEFAULTS.get(key, "NA"))
        try:
            config_text = build_params_config(values)
        except ValueError as exc:
            error = str(exc)
            config_text = ""

    return render_template("index.html", values=values, config_text=config_text, error=error)


@app.route("/download", methods=["POST"])
def download():
    values = DEFAULTS.copy()
    for key in DEFAULTS:
        if key == "addleadSNPs":
            values[key] = "1" if request.form.get("addleadSNPs") == "1" else "0"
        elif key in {"posMapAnnot", "xqtlsMapdss"}:
            values[key] = ":".join(request.form.getlist(key))
        else:
            values[key] = request.form.get(key, DEFAULTS.get(key, "NA"))
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
