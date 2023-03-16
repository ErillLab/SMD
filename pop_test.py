import timeit

if __name__ == "__main__":

    my_list = [45, 23, 18, 23, 8, 39, 23, 29, 2, 11, 78, 7,
               82, 33, 18, 28, 5, 54, 21, 61, 893, 20, 3]
    print(my_list)

    start = timeit.default_timer()

    for i in range(len(my_list)):
        j = i+1
        while j < len(my_list):
            if (my_list[i] - my_list[j]) < 8 and (my_list[i] - my_list[j]) > -8:
                print(my_list[j], j)
                my_list.remove(my_list[j])
                j -= 1
            j += 1

    stop = timeit.default_timer()

    print('Time: ', stop - start)

    print(my_list)

    my_list = [20, 12, 43, 6, 78, 30, 14, 23, 9, 30, 29, 10, 40, 62, 73, 4, 2, 35, 14, 22, 0, 10, 21, 80, 90,
               79, 54, 2, 34, 26, 56, 87, 3, 20, 65, 34, 1, 98, 75, 57, 45, 23, 18, 23, 8, 39, 23, 29, 2, 11, 78, 7,
               82, 33, 18, 28, 5, 54, 21, 61]
    print(my_list)

    '''
    start = timeit.default_timer()

    for i in range(len(my_list)):
        j = i+1
        while j < len(my_list):
            if (my_list[i] - my_list[j]) < 8 and (my_list[i] - my_list[j]) > -8:
                my_list[j] =
            j += 1
    '''

