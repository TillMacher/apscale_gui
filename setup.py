import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="apscale_gui",
    version="1.0.8",
    author="Till-Hendrik Macher",
    author_email="till-hendrik.macher@uni-due.de",
    description="Advanced Pipeline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data - Graphical User Interface",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/TillMacher/apscale_gui",
    packages=setuptools.find_packages(),
    license = 'MIT',
    install_requires = ['pySimpleGUI >= 4.15.2',
                        'apscale >= 1.4.2',
                        'demultiplexer >= 1.1.0',
                        'lastversion >= 1.3.3',
                        'boldigger >= 1.2.5',
                        'plotly >= 4.9.0',
                        'kaleido >= 0.2.1',
                        'xlsxwriter >= 3.0.2'
                        ],
    include_package_data = True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    entry_points = {
        "console_scripts" : [
            "apscale_gui = apscale_gui.__main__:main",
        ]
    },
)
