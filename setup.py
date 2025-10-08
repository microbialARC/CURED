from setuptools import setup, find_packages

setup(
    name="cured",
    version="1.0.0",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            "cured = cured.__init__:main",
        ],
    },
    python_requires=">=3.8",
)
