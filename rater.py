tweets = ['a','b','c']
scores = ['c' for i in range(len(tweets))]

if __name__ == '__main__':
    i = 0
    while i < len(tweets):
        scores[i] = input(tweets[i]+":\n")
        if scores[i] == 'b':
            i -= 1
        else:
            i += 1

    print(scores)
    
