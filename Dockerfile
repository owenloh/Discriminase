# Static hosting for the Discriminase web app (no backend).
# Caddy serves web/ with correct module MIME types + gzip, on Railway's $PORT.
FROM caddy:2-alpine
COPY web/ /srv/
COPY Caddyfile /etc/caddy/Caddyfile
# Caddy's default entrypoint runs /etc/caddy/Caddyfile; $PORT is provided by Railway.
