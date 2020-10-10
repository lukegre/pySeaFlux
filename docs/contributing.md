# Instructions for developers

## Development environment
I've tried to make the entry level to becoming a developer for SeaFlux 
as simple as possible. Here's a "tl;dr" version of what to do.

```bash
git clone https://github.com/luke-gregor/SeaFlux.git
cd SeaFlux

conda-env create --file ./dev-environment.yml
conda activate seaflux-dev

pre-commit install -t pre-commit
pre-commit install -t pre-push
```

This should be it for the setup.  
Here's a more indepth description of what the steps do:
- **Cloning** creates a local copy of the repo on your computer. When you `push` the results, it should automatically push to the SeaFlux repository (i.e. the upstream target is set to the correct URL).
- We then **create the environment** with conda from a file that specifies the required packages for development. There after that environment is activated. 
- Lastly we install **pre-commit** which will run syntax checks and tests before you record the changes made locally (commit) and send (push) your code to the original repository. 

## Creating a pull request for an feature/issue
I use **VS Code** for most of my git interaction. Their *Source Control* extension in combination with *Git Graph* and *Git Lense* extensions makes life a lot easier! I highly recommend this for those who do not want to spend hours watching Youtube tutorials on what `git rebase` means (I've been there). 

```bash
git checkout -b <issue_branch_name>  # creates new branch
```
Then make changes to the files you need to fix or add new features. 
```bash
# stashes changes for commit
git add <path_to_the_changed_file>  
# git commit will trigger syntax checks - you may have to add changes a second time
git commit -m "useful message"  # will store stashed changes with your comment
# git push will trigger tests
git push  # send the changes to GitHub
```

Lastly, Go to https://github.com/luke-gregor/SeaFlux and create a pull request


## Documentation
`mkdocs` is a very simple documentation package that uses markdown. It is quick and efficient, though not as powerful as `rst`.

* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs gh-deploy` - Build docs to gh-pages branch
* `mkdocs -h` - Print help message and exit.


## Deployment to PyPi
This is handled automatically when new tags are made. 
The version number of the package will be based on the tag name. 
