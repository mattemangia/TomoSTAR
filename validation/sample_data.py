# Copyright 2026 Matteo Mangiagalli
# SPDX-License-Identifier: Apache-2.0
"""
Real data for the validation of the readers, the filters, the responses and the picker: a few events
of the 2016-2017 central Italy sequence from the INGV FDSN web services (QuakeML with the bulletin
picks, StationXML with responses, miniSEED of the IV stations nearby). Downloaded once into the work
folder; the Norcia example's data are used instead when TOMOSTAR_NORCIA_DATA points to them.
"""
import os
import subprocess

from common import work

BASE = "https://webservices.ingv.it/fdsnws"
# Three events of different size: the Mw 6.5 Norcia mainshock, an Mw 4.x and an M 3.x aftershock.
EVENTS = {"8863681": "2016-10-30T06:40:17", "10403221": "2017-01-18T10:25:40", "10027911": "2016-11-16T11:52:20"}
STATIONS = ["NRCA", "CAMP", "MMO1", "FEMA", "ARRO", "CESI", "SSM1", "MNTP"]


def curl(url, path, data=None):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    cmd = ["curl", "-sS", "--compressed", "--retry", "4", "-o", path, url]
    if data:
        cmd[1:1] = ["--data-binary", data]
    subprocess.run(cmd, check=True)
    return path


def sample():
    """Paths of the QuakeML folder, the StationXML file and the miniSEED files."""
    norcia = os.environ.get("TOMOSTAR_NORCIA_DATA")
    if norcia:
        mseed = [os.path.join(norcia, "waveforms", e, f"{e}.mseed") for e in EVENTS if e != "8863681"]
        return os.path.join(norcia, "quakeml"), os.path.join(norcia, "stations.xml"), [m for m in mseed if os.path.exists(m)]
    folder = work("sample")
    qml = work("sample", "quakeml")
    for e in EVENTS:
        curl(f"{BASE}/event/1/query?eventid={e}&includearrivals=true&format=xml", os.path.join(qml, f"{e}.xml"))
    sta = curl(f"{BASE}/station/1/query?network=IV&station={','.join(STATIONS)}&starttime=2016-08-01&endtime=2017-03-01&level=response",
               os.path.join(folder, "stations.xml"))
    mseed = []
    for e, t in EVENTS.items():
        if e == "8863681":
            continue  # the mainshock clips the nearby records; the aftershocks are used for the waveforms
        path = os.path.join(folder, f"{e}.mseed")
        if not os.path.exists(path):
            import datetime as dt
            t0 = dt.datetime.fromisoformat(t)
            req = "\n".join(f"IV {s} * HH?,EH? {(t0 - dt.timedelta(seconds=20)).isoformat()} {(t0 + dt.timedelta(seconds=60)).isoformat()}" for s in STATIONS)
            open(path + ".req", "w").write(req + "\n")
            curl(f"{BASE}/dataselect/1/query", path, data="@" + path + ".req")
        mseed.append(path)
    return qml, sta, mseed
