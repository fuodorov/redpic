from setuptools import setup, find_packages
import redpic

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='redpic',
    version=redpic.__version__,
    author='Vyacheslav Fedorov',
    author_email='fuodorov1998@gmail.com',
    description=redpic.__doc__,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/fuodorov/redpic',
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics'
    ],
    python_requires='>=3.6',
)
