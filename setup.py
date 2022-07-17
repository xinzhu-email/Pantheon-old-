import setuptools

with open("README.md","r",encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "Patheon",
    version = "0.0.1",
    author = "Xinzhu Jiang",
    description = "A graphical interface of single cell analysis",
    long_description=long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/xinzhu-email/Pantheon",
    packages = setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Pyhton :: 3",
    ],
    python_requires = '>=3.8',
    install_requires = [
        "bokeh",
        "anndata",
        "pandas",
        "colorcet"
    ],
)