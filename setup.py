import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="apscale_gui",
    version="3.0.1",
    author="Till-Hendrik Macher",
    author_email="macher@uni-trier.de",
    description="Advanced Pipeline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data - Graphical User Interface",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/TillMacher/apscale_gui",
    packages=setuptools.find_packages(),
    license = 'MIT',
    install_requires = [
                        "apscale>=4.1.4",
                        "apscale_blast>=1.2.7",
                        "boldigger3>=2.1.4",
                        "demultiplexer2>=1.1.6",
                        "ete3>=3.1.3",
                        "lxml_html_clean>=0.4.2",
                        "playwright>=1.54.0",
                        "pandas>=2.3.1",
                        "numpy>=2.3.2",
                        "streamlit>=1.48.1",
                        "streamlit-file-browser>=3.2.22",
                        "scipy>=1.16.1",
                        "cutadapt>=5.1",
                        "update-checker>=0.18.0",
                        "powerlaw>=1.5",
                        "requests>=2.32.3",
                        "beautifulsoup4>=4.13.4",
                        "lxml>=6.0.0"
                        ],
    include_package_data = True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.12',
    entry_points = {
        "console_scripts" : [
            "apscale_gui = apscale_gui.__main__:main",
        ]
    },
)
