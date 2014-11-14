# a demonstration of buble sort and merge sort, p is a list of integers
# swap the values if p[i+1] < p[i]


def bubblesort(p):
    for j in range(len(p)):
        for i in range(len(p)-1):
           if p[i] > p[i+1]:
                 temp = p[i]
                 p[i] = p[i+1]
                 p[i+1] = temp
    return p

    # we can do better, when there is no swap, we stop the for loop


def bubblesort2(p):
      swapped = True
      while swapped:
            swapped = False
            for i in range(len(p)-1):
                  if p[i] > p[i+1]:
                        temp = p[i]
                        p[i] = p[i+1]
                        p[i+1] = temp
                        swapped = True
      return p


# divide and conqure method merge-sort

# divide the list in half, continue until we have single list, merge sub-lists

def mergesort(p):
      # print p
      if len(p) < 2 :
            return p[:]
      else:
            middle = len(p) / 2
            left = mergesort(p[:middle])   # this is the beauty of recurrsion, further break down the list to smaller lists
            right =mergesort(p[middle:])
            together = merge(left, right)
            #print "merged", together
            return together


      # now we need to define the merge function

def merge(left, right):
      result = []
      i,j = 0,0
      while i < len(left) and j < len(right):
                    if left[i] <= right[j]:
                        result.append(left[i])
                        i = i+1
                    else:
                        result.append(right[j])
                        j = j+1
      while i< len(left):
            result.append(left[i])
            i = i + 1
      while j < len(right):
            result.append(right[j])
            j = j + 1
      return result
