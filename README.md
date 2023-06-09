# the-friendly-stars
`the-friendly-stars` is a Python toolkit for interacting with catalogs of stars and images of the sky. It can be used to make custom finder charts, and it may soon include additional visualizations.

Following [The Friendly Stars](https://books.google.com/books?id=uUYJAAAAIAAJ&printsec=frontcover&source=gbs_ge_summary_r&cad=0#v=onepage&q&f=false) by Martha Evans Martin (1925), "the chief aim of this [package] is to share with others the pleasure which the writer has had in what may be called a relation of personal friendship with the stars." It is still very much so a *work in progress*.

*(If there's any chance you had been using a pre-2023 version of `the-friendly-stars`, please note it has moved to [`the-earlier-version-of-the-friendly-stars`](https://github.com/zkbt/the-earlier-version-of-the-friendly-stars). We're starting from scratch here!)*

## Installation
You should be able to install this simply by running `pip install git+https://github.com/zkbt/the-friendly-stars.git`.

If you want to be able to modify the code yourself, please also feel free to fork/clone this repository onto your own computer and install directly from that editable package. For example, this might look like:
```
git clone git@github.com:zkbt/the-friendly-stars.git
cd the-friendly-stars/
pip install -e '.[develop]'
```
This will link the installed version of the `thefriendlystars` package to your local repository. Changes you make to the code in the repository should be reflected in the version Python sees when it tries to `import thefriendlystars`. The `-e` means "we can edit the code in its local director"y; the `.` means "install from this current directory"; and the `[develop]` means "install the development dependencies for testing and documentation".

## Contributors

Code or conceptual contributions were made by:
- [Zach Berta-Thompson](https://github.com/zkbt)
- [Luci Ibarra Perez](https://github.com/luib0557)
- (maybe you soon!)

The [ESA Gaia mission](https://gea.esac.esa.int/archive/) provides the data that sits at the core of `the-friendly-stars`. Any work that makes use of this tiny package should acknowledge the incredible work of the thousands of people who make Gaia possible and cite the Gaia archive.
