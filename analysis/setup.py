from setuptools import setup

setup(
    name="bbackend",
    version="0.0.0",
    author="Yan-Rong Li",
    packages=["bbackend",],
    package_dir={'bbackend':''},
    install_requires=["numpy","matplotlib","scipy"],
)