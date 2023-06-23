from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="supernova-cffi",
    version="1.0.0",
    author="Mingde Yin",
    description="Orbit propagation and analysis tool.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/spacesys-finch/supernova/",
    license="LGPLv3",
    packages=["supernova"],
    setup_requires=["cffi>=1.15.0"],
    cffi_modules=["supernova/_builder.py:ffibuilder"],
    install_requires=["cffi>=1.15.0", "numpy", "matplotlib"],
    extras_require={"test": ["pytest", "poliastro"], "demo": ["celest"]},
    options={"build_ext": {"inplace": True}},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: C",
        "Topic :: Scientific/Engineering :: Physics",
        "Intended Audience :: Science/Research",
    ],
    python_requires=">=3.8",
)
