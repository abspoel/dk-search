import argparse
import uuid
import subprocess

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bits", type=int)
    parser.add_argument("d0", type=int)
    parser.add_argument("d1", type=int)
    parser.add_argument("num_parts", type=int)
    parser.add_argument("start", type=int)
    parser.add_argument("stop", type=int)
    parser.add_argument("host")
    parser.add_argument("-k", type=int, default=1)
    parser.add_argument("--no-star", action='store_true')
    parser.add_argument("--cwd", default="dk-search")
    parser.add_argument("--cmd", default="./target/release/dk-search")
    parser.add_argument("--no-ssh", action='store_true')
    parser.add_argument("--log-dir", default="output")

    args = parser.parse_args()

    for part in range(args.start, args.stop):
        print("Part {} (total={} parts)".format(part, args.num_parts))
        max_p = (1 << args.bits) - 1
        start = round(part * (max_p / args.num_parts))
        stop = round((part + 1) * (max_p / args.num_parts))
        if stop > max_p:
            stop = max_p
        star = "--dk-star" if not args.no_star else ""
        params = "-k {} {} --start {} {} {} {}".format(args.k, star, start, args.d0, args.d1, stop)
        desig = "{}-{}{}{}".format(args.d0, args.d1, "-d2" if args.k == 2 else "", "-star" if not args.no_star else "")
        filename = "{bits}bits-{desig}-{num_parts}-parts--{part}".format(bits=args.bits, desig=desig, num_parts=args.num_parts, part=part)
        command = "nohup {cmd} {params} >{log_dir}/{filename}.log 2>{log_dir}/{filename}.err </dev/null &".format(cwd=args.cwd, cmd=args.cmd, log_dir=args.log_dir, params=params, filename=filename)
        if not args.no_ssh:
            command = ["ssh", args.host, "cd {cwd}; " + command]
            shell = False
        else:
            shell = True
        print(command)
        subprocess.run(command, shell=shell)

if __name__ == "__main__":
    main()
