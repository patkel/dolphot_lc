# How to Contribute
Thank you for considering contributing to Dolphot-LC! You can contribute in several ways.

## Opening issues

If you find a bug or wish to propose a new feature, open a new issue on GitHub. Clearly write your proposal for change and include all relevant information. 

## Discussing the project

If you have questions on how to use the package or share relevant research, open a Discussion on GitHub. Feel free to also answer in other users discussions. 

## Contributing software

All Dolphot-LC code is hosted on GitHub. You can improve this package by proposing new functionality, solving bugs or implementing accepted feature requests.
Please ensure you own the rights to your code and open a discussion first to prevent duplicate efforts. 


### Make a fork (copy) of Dolphot_LC

**You only need to do this once**

1. Go to the [Dolphot_LC repository home page](https://github.com/patkel/dolphot_lc)
3. Click on the *Fork* button (top-right-hand corner)
4. Select the namespace that you want to create the fork in, this will usually be your personal namespace

### Clone your fork

```bash
git clone https://github.com/<username>/dolphot_lc.git
```

### Updating your fork

If you already have a fork of Dolphot_LC, and are starting work on a new project you can link your clone to the main (`patkel`) repository and pull in changes that have been merged since the time you created your fork, or last updated:

1. Link your fork to the main repository:

    ```bash
    cd dolphot_lc
    git remote add patkel https://github.com/patkel/dolphot_lc.git
    ```
 

2. Fetch new changes from the `patkel` repo

    ```bash
    git fetch patkel
    ```

### Creating a new feature branch

All changes should be developed on a feature branch, in order to keep them separate from other work, simplifying review and merge once the work is done.

To create a new feature branch:

```bash
git fetch patkel
git checkout -b my-new-feature patkel/master
```

### Hack away

1. Develop the changes you would like to introduce, using `git commit` to finalise a specific change.
   Ideally commit small units of change often, rather than creating one large commit at the end, this will simplify review and make modifying any changes easier.

    Commit messages should be clear, identifying which code was changed, and why.
   Common practice is to use a short summary line (<50 characters), followed by a blank line, then more information in longer lines.

2. Push your changes to the remote copy of your fork on GitHub

    ```bash
    git push origin my-new-feature
    ```
   **Note:** For the first `push` of any new feature branch, you will likely have to use the `-u/--set-upstream` option to `push` to create a link between your new branch and the `origin` remote:

    ```bash
    git push --set-upstream origin my-new-feature
    ```

### Open a Pull Request

When you feel that your work is finished, you should create a Pull Request to propose that your changes be merged into the main (`patkel`) repository.

After you have pushed your new feature branch to `origin`, you should find a new button on the [Dolphot_LC repository home page](https://github.com/patkel/dolphot_lc) inviting you to create a Pull Request out of your newly pushed branch.
You should click the button, and proceed to fill in the title and description boxes on the PR page.

Once the request has been opened, one of the maintainers will review the change.

