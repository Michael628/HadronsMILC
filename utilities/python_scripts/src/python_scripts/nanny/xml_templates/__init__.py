def wrapper(runid: str, sched: str, cfg: str) -> dict:
    params = {
        "grid": {
            "parameters": {
                "runId": runid,
                "trajCounter": {
                    "start": cfg,
                    "end": "10000",
                    "step": "10000",
                },
                "genetic": {
                    "popSize": "20",
                    "maxGen": "1000",
                    "maxCstGen": "100",
                    "mutationRate": "0.1",
                },
                "graphFile": "",
                "scheduleFile": sched,
                "saveSchedule": "false",
                "parallelWriteMaxRetry": "-1",
            },
            "modules": {},
        },
    }

    return params
