import threading as thg

def rev(thread):
    thread.join()
    print('Поток завершил работу.')
    return 1