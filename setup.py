import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="metaprocessor",
    version="0.3.beta",
    author="Till-Hendrik Macher",
    author_email="till-hendrik.macher@uni-due.de",
    description="Metaprocessor - a GUI-based, platform-independent (e)DNA metabarcoding processing pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/metaprocessor",
    packages=setuptools.find_packages(),
    license = 'MIT',
    install_requires = ['pySimpleGUI>=4.15.2',
                        'pandas>=0.25.3',
                        'numpy>=1.18.1',
                        'xlrd>=1.2.0',
                        'openpyxl>=3.0.3',
                        'xlsxwriter>=1.2.7',
                        'biopython>=1.77',
                        'plotly>=4.9.0',
                        'kaleido>=0.0.3',
                        'cutadapt>=2.1'
                        'openpyxl >= 3.0.3',
                        'psutil >= 5.7.3',
                        'joblib >= 0.16.0',
                        ],
    include_package_data = True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
