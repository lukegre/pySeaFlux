from setuptools import setup


setup(
    include_package_data=True,
    use_scm_version={
        "write_to": "seaflux/_version.py",
        "write_to_template": '__version__ = "{version}"',
        "tag_regex": r"^(?P<prefix>v)?(?P<version>[^\+]+)(?P<suffix>.*)?$",
    },
)
