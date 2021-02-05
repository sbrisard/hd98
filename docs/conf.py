import pathlib

project_root = pathlib.Path.cwd().parent


def get_metadata(key):
    with open(project_root / "metadata" / (key + ".txt"), "r", encoding="utf8") as f:
        return f.read()


project = "hd98"
copyright = get_metadata("year") + ", " + get_metadata("author")
author = get_metadata("author")
release = get_metadata("version")


extensions = ["sphinx.ext.autodoc", "sphinx.ext.todo", "breathe"]

templates_path = ["_templates"]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

todo_include_todos = True


breathe_projects_source = {
    "hd98": (
        str(project_root / "include" / "hd98"),
        ["hd98.hpp"],
    )
}

breathe_doxygen_config_options = {"GENERATE_TODOLIST": "YES"}

html_theme = "alabaster"
html_static_path = ["_static"]
