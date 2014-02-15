from kombu import Queue

BROKER_URL = 'amqp://guest@godel97'
CELERY_RESULT_BACKEND = 'redis://godel97/0'
CELERY_TASK_SERIALIZER = 'pickle'
CELERY_RESULT_SERIALIZER = 'pickle'
CELERY_ACCEPT_CONTENT=['pickle', 'json']
CELERY_TIMEZONE = 'America/New York'
CELERY_ENABLE_UTC = True
CELERYD_PREFETCH_MULTIPLIER = 100

CELERY_QUEUES = (
    Queue('celery', routing_key='celery'),
    Queue('transient', routing_key='transient', delivery_mode=1))
