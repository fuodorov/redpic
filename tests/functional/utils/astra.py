import pandas as pd

from redpic import constants as const


def read_beam_astra_particles(file_name):
    cols = [
        "x",
        "y",
        "z",
        "px",
        "py",
        "pz",
        "clock",
        "charge",
        "id",
        "flag",
    ]  # m    m    m    eV/c  eV/c  eV/c  ns       nC
    df = pd.read_csv(file_name, header=None, delim_whitespace=True, names=cols, dtype="float32")
    df = df[df.flag != -15]  # ignore the lost particles
    df["px"] = df["px"] / 1e6  # MeV/c
    df["py"] = df["py"] / 1e6  # MeV/c
    df0 = df.head(1)  # remove the reference particle
    df = df.drop(df0.index)
    z0 = df0.z.values[0]
    pz0 = df0.pz.values[0]
    df["z"] = z0 + df["clock"] * 1e-9 * const.c  # m
    df["pz"] = (pz0 + df["pz"]) / 1e6  # MeV/c
    return df


def read_track_astra_particles(file_name):
    cols = ["x", "y", "z", "px", "py", "pz", "clock", "charge", "id", "flag"]
    df = pd.read_csv(file_name, header=None, delim_whitespace=True, names=cols, dtype="float32")
    df = df[df.flag != -15]  # ignore the lost particles
    df["px"] = df["px"] / 1e6  # MeV/c
    df["py"] = df["py"] / 1e6  # MeV/c
    df0 = df.head(1)
    df = df.drop(df0.index)
    z0 = df0.z.values[0]
    pz0 = df0.pz.values[0]
    df["z"] = z0 + df["z"]  # m
    df["pz"] = (pz0 + df["pz"]) / 1e6  # MeV/c
    return df
