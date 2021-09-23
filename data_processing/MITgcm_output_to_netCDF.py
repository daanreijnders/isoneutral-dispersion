import xarray as xr
import xmitgcm
import argparse
import os
import sys


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Convert MITgcm native output to netCDF using `xmitgcm`."
    )
    parser.add_argument("path", type=dir_path, help="Path with .mds files.")
    parser.add_argument("-n", "--name", help="Name of output file", type=str)
    parser.add_argument("-t0", help="Initial timestamp (`YYYY-MM-DD`)", default="1900-01-01", type=str)
    parser.add_argument("-dt", help="Number of seconds per timestep", default=300, type=int)
    parser.add_argument(
        "--addSalinity",
        action="store_true",
        help="Add uniform (35 psu) salinity field.",
    )
    return parser.parse_args()


def dir_path(path):
    if os.path.isdir(path):
        return os.path.realpath(path)
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == "":
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' " "(or 'y' or 'n').\n")


def validate(parsed_args):
    assert len(parsed_args.t0) == 10, "`t0` must be formatted as `YYYY-MM-DD`."
    assert parsed_args.t0[4] == '-' and parsed_args.t0[7]== '-' , "`t0` must be formatted as `YYYY-MM-DD`."
            
    
def main():
    parsed_args = parse_arguments()
    print("Input path:", parsed_args.path)
    print("Output path:", parsed_args.name)
    print("Adding uniform salinity:", parsed_args.addSalinity)
    
    validate(parsed_args)

    if parsed_args.name:
        outpath = parsed_args.path + "/" + parsed_args.name
    else:
        outpath = parsed_args.path + "/" + "out.nc"
    if outpath[-3:] != ".nc":
        outpath = outpath + ".nc"

    write = True
    if os.path.exists(outpath):
        if query_yes_no("File already exists. Overwrite?"):
            pass
        else:
            write = False
    if write:
        ds = xmitgcm.open_mdsdataset(
            parsed_args.path,
            prefix=["surfDiag", "dynDiag", "diffElem"],
            read_grid=True,
            geometry="cartesian",
        )
        # Correct timestamps
        ds = ds.assign_coords(
            {"time": ds["time"].values.astype("timedelta64[s]").astype("float64") * parsed_args.dt}
        )
        
        t0 = parsed_args.t0
        t0_y, t0_m, t0_d = t0[0:4], t0[5:7], t0[8:]
        
        ds["time"].attrs["units"] = f"seconds since {t0_y}-{t0_m}-{t0_d}"
        ds["time"].attrs["calendar"] = "360_day"
        ds = xr.decode_cf(ds)

        if parsed_args.addSalinity:
            SALT = (
                (ds.THETA * 0 + 35.0)
                .compute()
                .rename("SALT")
                .assign_attrs(
                    {"standard_name": "SALT", "long_name": "Salinity", "units": "psu"}
                )
            )
            ds["SALT"] = SALT
        ds.to_netcdf(outpath)
        
    else:
        pass


if __name__ == "__main__":
    main()
