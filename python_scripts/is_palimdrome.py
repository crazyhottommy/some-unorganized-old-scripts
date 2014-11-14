#! /usr/bin/env python 


#This script demonstrate recursion
#check if a given string is a palindrome


def is_palindrome(s):
    ''' (str) ----> bool
    Return Ture or False if a str is palindrome or not
    >>> is_palindrome("a")
    True
    >>> is_palindrome("ab")
    False
    >>> is_palindrome("abba")
    True
    '''
    if len(s) <= 1: 
        return True
    else:
        return s[0]==s[-1] and is_palindrome(s[1:-1])
import doctest
doctest.testmod()

#Fibonacci numbers
#0  1 1 2 3 5 8 13 21....

def fib(x):
    '''Retrun fibonacci of x, where x is a non-negative int'''
    if x==0 or x==1:
        return 1
    else:
        return fib(x-1) + fib(x-2)
    
# pythonic way 

def is_palindromev2(s):
    return s==s[::-1]
