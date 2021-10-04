import os

# gunicorn config file
# see :
# http://docs.gunicorn.org/en/stable/configure.html
# http://docs.gunicorn.org/en/stable/settings.html

bind = "0.0.0.0:{}".format(int(os.getenv("PORT", 8000)))
workers = 1
threads = 4
errorlog = "logs/gunicorn-error.log"
accesslog = "logs/gunicorn-access.log"
loglevel = os.getenv("GUNICORN_LOG_LEVEL", "info")
