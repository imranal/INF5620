def fib(n):
    """Print the Fibonacci series up to n."""
    a, b = 0, 1
    while b < n:
        print "%g"%b,
        a, b = b, a + b

if __name__ == "__main__":
    fib(8139748913419874189489179471297491718749812980000041890174124104812001274012401927000000000000)
