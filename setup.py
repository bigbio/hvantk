from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="pyvatk",
    version="0.1.0",
    description="A package for gene and variant annotation.",
    author="Enrique Audain & Rafiga Masmaliyeva",
    author_email="enrique.audain.martinez@uni-oldenburg.de",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",  # TODO: Update github url
    packages=find_packages(),
    install_requires=[
        "click",
        "setuptools",
        "hail"
    ],
    entry_points={
        'console_scripts': [
            'pyvatk=pyvatk.pyvatk_cli:pyvatk_main',
        ],
    },
    classifiers=[
        # Classifiers for your package
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires='>=3.10.0',
)