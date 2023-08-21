from setuptools import setup
from pathlib import Path

TOP_DIR = Path(__file__).parent.resolve()

with open(TOP_DIR.joinpath("requirements.txt")) as handle:
    REQUIREMENTS = [f"{_.strip()}" for _ in handle.readlines() if " " not in _]

setup(name='environmentfinder',
      version='1.1',
      description='Tool for finding atomic environments in crystal structures',
      url='http://github.com/PabloPiaggi/EnvironmentFinder',
      author='Pablo Piaggi',
      author_email='ppiaggi@princeton.edu',
      license='GNU General Public License v3.0',
      packages=['environmentfinder'],
      install_requires=REQUIREMENTS,
      scripts=['bin/environmentfinder'],
      include_package_data=True,
)
