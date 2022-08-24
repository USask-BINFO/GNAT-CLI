# AcrVis

Generalizable tool used to visualize the genomic neighbourhoods of Acr/Aca proteins and their taxonomic data.

## Installation

To install the neccessary scripts needed to run the pipeline, simply run the following command in the directory you wish to have the project installed:

```
git clone https://git.cs.usask.ca/njt694/acrvis.git
```

If you wish to install any updates to the pipeline, cd to the directory created by the git clone and run:

```
git pull
```

Then add the scripts/ folder to your PATH, if you don't want to specify the location every time you call main.sh. 
You can add the scripts to your PATH on unix systems with the following command (assuming you git cloned into your home directory):

```
export PATH="~/acrvis/scripts:$PATH"
```

## Database

To run the pipeline, a database is required for psi-BLAST searches and building the genomic neighbourhoods. Scripts to generate/update your own databases are located in the db_scripts/ folder.

To build your own database you must specify a link to the ncbi ftp site. The link should start with "ftp://ftp.ncbi.nih.gov/genomes/genbank/" and should end at an "assembly_summary.txt" file. For example, to build a database of all the bacterial genomes found on ncbi in your current directory, you would run:

```
db_scripts/get_ncbi_genomes.sh -l "ftp://ftp.ncbi.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"
```

To update that database, assuming it was created in your home directory, you would run:

```
db_scripts/get_ncbi_genomes.sh -l "ftp://ftp.ncbi.nih.gov/genomes/genbank/bacteria/assembly_summary.txt" -e ~/bacterial_db
```

For information on all of the optional arguments when building a database, run:

```
db_scripts/get_ncbi_genomes.sh -h
```

Note that some of these databases may be very large (terabytes!), so make sure the machine you are downloading them onto has enough space. Also due to the large size, the script may take a long time to run, so it is recommended to run it inside of a 'screen'. You can start and attach to a new screen with the name 'db_building' with the following code:

```
screen -S db_building
```

To detach from the screen, use ctrl+a, and to re-attach to the screen, run:

```
screen -r db_building
```

## Usage


## Name
Choose a self-explaining name for your project.

## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might be unfamiliar with. A list of Features or a Background subsection can also be added here. If there are alternatives to your project, this is a good place to list differentiating factors.

## Badges
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.

## Visuals
Depending on what you are making, it can be a good idea to include screenshots or even a video (you'll frequently see GIFs rather than actual videos). Tools like ttygif can help, but check out Asciinema for a more sophisticated method.

## Installation
Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew. However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a specific context like a particular programming language version or operating system or has dependencies that have to be installed manually, also add a Requirements subsection.

## Usage
Use examples liberally, and show the expected output if you can. It's helpful to have inline the smallest example of usage that you can demonstrate, while providing links to more sophisticated examples if they are too long to reasonably include in the README.

## Support
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.

## Roadmap
If you have ideas for releases in the future, it is a good idea to list them in the README.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.

## Temp Notes
If getting AttributeError: 'NoneType' object has no attribute 'split' when running find_clusters.py, mentioning 'threadpoolctl.py', try running pip install threadpoolctl.py -U to update the package. The error doesn't prevent the pipeline from running correctly, but clutters output.
