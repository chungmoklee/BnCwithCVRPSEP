from setuptools import setup, find_namespace_packages

setup(
    name = 'cvrpsep',
    version = '0.1',    
    setup_requires=["cffi>=1.0.0"],
    packages = find_namespace_packages(where='.'),
    cffi_modules = ["cvrpsep/buildcvrpsep.py:ffibuilder"],
    install_requires = ["cffi>=1.0.0"],
    package_data = {'cvrpsep': ['Cvrpsep.pdf'], 'cvrpsep.src': ['*']}
)