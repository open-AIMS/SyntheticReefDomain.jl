# Registers Pluto with jupyter-server-proxy on Binder (copied to ~/.jupyter by postBuild),
# so it is reachable at <binder url>/pluto/ and from the JupyterLab launcher.
c.ServerProxy.servers = {
    "pluto": {
        "command": ["julia", "--project=notebooks", "binder/start_pluto.jl", "{port}"],
        "timeout": 300,
        "launcher_entry": {"title": "Pluto.jl (reef builder)"},
    }
}
