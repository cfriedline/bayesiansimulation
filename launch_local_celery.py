from fabric.operations import local
import multiprocessing as mp
import time

rabbit = "/Users/chris/src/rabbitmq_server-3.2.3/sbin/rabbitmq-server"
worker = "celery worker -A runsimulation4 -l info --autoscale=1,1"


def run_task(task):
    return local(task)


def main():
    jobs = []
    pool = mp.Pool(2)
    jobs.append(pool.apply_async(run_task, args=(rabbit,)))
    time.sleep(5)
    jobs.append(pool.apply_async(run_task, args=(worker,)))

    while (True):
        for j in jobs:
            print j.get()
        time.sleep(2)


if __name__ == '__main__':
    main()