import argparse
import os
import shutil
import sys
from glob import glob

file_patterns = ["datacard_ggH_4g_m*_ct*_*_*.txt","datacardInputs_ggH_4g_m*_ct*_*_*.root"]
directory_patterns = ["m*_ct*_*_*_4g_ggH"]
def dir_size(path):
    total=0
    for root, _, files in os.walk(path):
        for f in files:
            fp = os.path.join(root, f)
            if not os.path.islink(fp):
                total += os.path.getsize(fp)
    return total
def find_targets(directory):
    files,dirs = [],[]
    for pattern in file_patterns:
        for path in glob(os.path.join(directory, pattern)):
            if os.path.isfile(path):
                files.append((path, os.path.getsize(path)))
    for pattern in directory_patterns:
        for path in glob(os.path.join(directory, pattern)):
            if os.path.isdir(path) and not os.path.islink(path):
                dirs.append((path, dir_size(path)))
    return sorted(set(files)), sorted(set(dirs))





def report(files, dirs, directory):
    teal="\033[38;5;44m"
    mustard="\033[38;5;136m"
    reset="\033[0m"
    print(f"{teal}Datacard outputs found in {directory}:{reset}")
    for label,items in (("file",files), ("directory",dirs)):
        if not items:
            continue
        print(f"{mustard}  {len(items)} {label}{'' if len(items) == 1 else 's'}:{reset}")
        for path,size in items:
            suffix="/" if label=="directory" else ""
            print(f"    {os.path.basename(path)}{suffix}  ")

def main():
    parser= argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", dest="directory", default=os.path.dirname(os.path.abspath(__file__)),help="directory to clean (default: the directory holding this script)")
    parser.add_argument("-n", "--dry-run", action="store_true",help="list what would be deleted and exit without deleting")
    parser.add_argument("-y", "--yes", action="store_true",help="delete without asking for confirmation")
    args = parser.parse_args()
    green = "\033[1;92m"
    red = "\033[1;91m"
    reset = "\033[0m"
    directory = os.path.abspath(args.directory)
    if not os.path.isdir(directory):
        sys.exit(f"{red}no such directory: {directory}{reset}")

    files,dirs=find_targets(directory)
    if not files and not dirs:
        print(f"{green}Nothing to clean in {directory}.{reset}")
        return
    report(files,dirs,directory)
    if args.dry_run:
        print(f"{green}Dry run -- nothing deleted.{reset}")
        return

    if not args.yes:
        answer = input(f"{red}Delete these permanently? [y/N] {reset}").strip().lower()
        if answer not in ("y", "yes"):
            print("Aborted, nothing deleted.")
            return
    removed, failed = 0,0
    for path, _ in files:
        try:
            os.remove(path)
            removed += 1
        except OSError as e:
            print(f"{red}  could not remove {path}: {e}{reset}")
            failed += 1
    for path, _ in dirs:
        try:
            shutil.rmtree(path)
            removed += 1
        except OSError as e:
            print(f"{red}  could not remove {path}: {e}{reset}")
            failed += 1
    print(f"{green}Removed {removed} item{'' if removed == 1 else 's'}")
    if failed:
        sys.exit(f"{red}{failed} items could not be removed{reset}")

if __name__ == "__main__":
    main()
