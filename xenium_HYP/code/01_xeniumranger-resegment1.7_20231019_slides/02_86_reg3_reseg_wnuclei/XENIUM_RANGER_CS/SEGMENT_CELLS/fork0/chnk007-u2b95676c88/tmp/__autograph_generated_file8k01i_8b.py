# coding=utf-8
from __future__ import annotations
def outer_factory():

    def inner_factory(ag__):

        def tf__printwhen(*args):
            """Print date and time along with message."""
            with ag__.FunctionScope('printwhen', 'fscope', ag__.STD) as fscope:
                now = ag__.converted_call(ag__.ld(datetime).datetime.now, (), None, fscope)
                pid = ag__.converted_call(ag__.ld(os).getpid, (), None, fscope)
                rusage = ag__.converted_call(ag__.ld(resource).getrusage, (ag__.ld(resource).RUSAGE_SELF,), None, fscope)
                rss_gb = ag__.ld(rusage).ru_maxrss / 1024 ** 2
                cpu_time = ag__.ld(rusage).ru_utime + ag__.ld(rusage).ru_stime
                cpu_factor = 1
                timestamp = ag__.converted_call(ag__.ld(now).timestamp, (), None, fscope)

                def get_state_2():
                    return (cpu_factor,)

                def set_state_2(vars_):
                    nonlocal cpu_factor
                    (cpu_factor,) = vars_

                def if_body_2():
                    nonlocal cpu_factor
                    start_timestamp = ag__.ld(timestamp) - ag__.ld(CPU_INTERVAL)
                    index = ag__.converted_call(ag__.ld(next), ((ag__.ld(i) for (i, x) in ag__.converted_call(ag__.ld(enumerate), (ag__.ld(CPU_TIMES),), None, fscope) if ag__.ld(x)[0] > ag__.ld(start_timestamp)), ag__.converted_call(ag__.ld(len), (ag__.ld(CPU_TIMES),), None, fscope)), None, fscope)

                    def get_state():
                        return ()

                    def set_state(block_vars):
                        pass

                    def if_body():
                        del ag__.ld(CPU_TIMES)[:ag__.ld(index) - 1]

                    def else_body():
                        pass
                    ag__.if_stmt(ag__.ld(index) > 1, if_body, else_body, get_state, set_state, (), 0)

                    def get_state_1():
                        return (cpu_factor,)

                    def set_state_1(vars_):
                        nonlocal cpu_factor
                        (cpu_factor,) = vars_

                    def if_body_1():
                        nonlocal cpu_factor
                        cpu_factor = (ag__.ld(cpu_time) - ag__.ld(CPU_TIMES)[0][1]) / (ag__.ld(timestamp) - ag__.ld(CPU_TIMES)[0][0])

                    def else_body_1():
                        nonlocal cpu_factor
                        pass
                    ag__.if_stmt(ag__.ld(timestamp) > ag__.ld(CPU_TIMES)[0][0], if_body_1, else_body_1, get_state_1, set_state_1, ('cpu_factor',), 1)

                def else_body_2():
                    nonlocal cpu_factor
                    pass
                index = ag__.Undefined('index')
                start_timestamp = ag__.Undefined('start_timestamp')
                ag__.if_stmt(ag__.ld(CPU_TIMES), if_body_2, else_body_2, get_state_2, set_state_2, ('cpu_factor',), 1)
                ag__.converted_call(ag__.ld(CPU_TIMES).append, ((ag__.ld(timestamp), ag__.ld(cpu_time)),), None, fscope)
                ag__.ld(print)(ag__.ld(now), f'[{ag__.ld(pid)}; {ag__.ld(rss_gb):.3f} GiB; {100 * ag__.ld(cpu_factor):.0f}%]', *ag__.ld(args))
                ag__.converted_call(ag__.ld(sys).stdout.flush, (), None, fscope)
        return tf__printwhen
    return inner_factory