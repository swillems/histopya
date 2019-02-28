#!venv/bin/python


import multiprocessing as mp


def parallelizedGenerator(
    function,
    function_args={},
    process_count=1,
):
    '''TODO comment'''
    out_queue = mp.Queue()
    processes = []
    for process_index in range(process_count):
        process_args = {
            "process_index": process_index,
            "process_count": process_count,
            "out_queue": out_queue,
        }
        process_args.update(function_args)
        process = mp.Process(
            target=function,
            args=(process_args, )
        )
        processes.append(process)
        process.start()
    running_processes = process_count
    while True:
        result = out_queue.get()
        if result is None:
            running_processes -= 1
            if running_processes == 0:
                break
        else:
            yield result
    while processes:
        process = processes.pop()
        process.join()
        process.terminate()


def partitionedQueue(data_to_partition, process_count, close_with_none=True):
    '''TODO comment'''
    queue = mp.Queue()
    for i in range(process_count):
        queue.put(range(i, len(data_to_partition), process_count))
    if close_with_none:
        queue.put(None)
    return queue


if __name__ == '__main__':
    pass














def serialGenerator(
    function,
    data_to_act_on,
    function_args={},
    process_count=1,
):
    '''TODO comment'''
    tally = len(data_to_act_on)
    out_queue = mp.Queue()
    in_queue = mp.Queue()
    # TODO in_queue.put()
    function_args["process_count"] = process_count
    function_args["out_queue"] = out_queue
    function_args["in_queue"] = in_queue
    for process_index in range(process_count):
        function_args["process_index"] = process_index
        process = mp.Process(
            target=function,
            args=(function_args, )
        )
        process.start()
    running_processes = process_count
    while running_processes > 0:
        to_reprocess, result = out_queue.get()
        if to_reprocess:
            in_queue.put(result)
        else:
            if result is None:
                running_processes -= 1
            else:
                yield result
                tally -= len(result)
                if tally == 0:
                    for process_index in range(process_count):
                        in_queue.put(None)


def parallelGenerator(
    function,
    data_to_act_on,
    function_args={},
    process_count=1,
):
    '''TODO comment'''
    out_queue = mp.Queue()
    function_args["out_queue"] = out_queue
    function_args["process_count"] = process_count
    for process_index in range(process_count):
        function_args["process_index"] = process_index
        function_args["data_range"] = range(process_index, len(data_to_act_on), process_count)
        process = mp.Process(
            target=function,
            args=(function_args, )
        )
        process.start()
    for process_index in range(process_count):
        yield out_queue.get()
