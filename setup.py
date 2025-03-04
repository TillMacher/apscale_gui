import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="apscale_gui",
    version="2.0.5",
    author="Till-Hendrik Macher",
    author_email="macher@uni-trier.de",
    description="Advanced Pipeline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data - Graphical User Interface",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/TillMacher/apscale_gui",
    packages=setuptools.find_packages(),
    license = 'MIT',
    install_requires = [
                        'apscale>=3.0.0',
                        'apscale_blast>=1.0.2',
                        'boldigger3>=1.1.1',
                        'pandas>=2.2.3',
                        'update_checker>=0.18.0',
                        'xlsxwriter>=3.2.0'
                        ],
    include_package_data = True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
    entry_points = {
        "console_scripts" : [
            "apscale_gui = apscale_gui.__main__:main",
        ]
    },
)
