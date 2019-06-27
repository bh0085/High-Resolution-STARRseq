import os, json, re
import pandas as pd
from _config import LOGS_PLACE

def get_last_runtime(NAME):
    logs_dir = os.path.join(LOGS_PLACE,NAME)
    files = os.listdir(logs_dir)
    unique_runtimes = [
        int(float(re.compile("^._([\d\.]+)").search(fn).group(1)))
        for fn in files
        if re.compile("^._[\d\.]+").search(fn)]
    recent_runtime = max(unique_runtimes)
    return recent_runtime

def read_last_output_logs(NAME):
    logs_dir = os.path.join(LOGS_PLACE,NAME)
    time = get_last_runtime(NAME)
    o_logs = [f for f in os.listdir(logs_dir) if f[0]=="o" and str(time) in f]

    #etexts = 
    outputs = pd.DataFrame()
    for fn in o_logs:
        with open(os.path.join(logs_dir,fn)) as fopen:
            txt = fopen.read()
            match =  re.compile("<json>(.*)</json>").search(txt)
            if match:
                json_format =match.group(1)
                jdata = json.loads(json_format)
                outputs = outputs.append(pd.Series(jdata).rename(fn))
            else: 
                outputs = outputs.append(pd.Series().rename(fn))
    return outputs

def read_last_error_logs(NAME):
    logs_dir = os.path.join(LOGS_PLACE,NAME)
    time = get_last_runtime(NAME)
    e_logs = [f for f in os.listdir(logs_dir) if f[0]=="e" and str(time) in f]
    output = []
    for fn in e_logs:
        with open(os.path.join(logs_dir,fn)) as fopen:
            output.append(fopen.read())
    return output