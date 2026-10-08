# Starts Pluto behind jupyter-server-proxy on Binder (see binder/jupyter_server_config.py).
#   julia --project=notebooks binder/start_pluto.jl <port>
import Pluto

port = parse(Int, ARGS[1])
Pluto.run(;
    host="0.0.0.0",
    port=port,
    launch_browser=false,
    # Binder already controls access to the session; Pluto's own secret would block the
    # proxied links (as in https://github.com/fonsp/pluto-on-binder)
    require_secret_for_access=false,
    require_secret_for_open_links=false,
    warn_about_untrusted_code=false,
    dismiss_update_notification=true
)
