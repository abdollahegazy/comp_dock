import subprocess
from pathlib import Path

species = ["DouglasFir","Arabidopsis","Eucalyptus","Human"]

mol_list = [p.stem for p in Path('../SmallMolecules/pdbfiles').glob('*.pdb')]


def run_command(command):
    print(f"Running: {command}")
    subprocess.run(command, shell=True, check=True,stdout=subprocess.DEVNULL)


def run_command_in_conda(command, env):
    full_command = (
        f"conda run -n {env} "
        f"{command}"
    )
    run_command(full_command)


