#!/usr/bin/env python3
"""Simple calculator supporting basic arithmetic operations."""


def add(a: float, b: float) -> float:
    return a + b


def subtract(a: float, b: float) -> float:
    return a - b


def multiply(a: float, b: float) -> float:
    return a * b


def divide(a: float, b: float) -> float:
    if b == 0:
        raise ValueError("Division by zero")
    return a / b


def main():
    ops = {"+": add, "-": subtract, "*": multiply, "/": divide}
    print("Simple Calculator (e.g. 2 + 3)")
    while True:
        try:
            line = input("> ").strip()
            if not line:
                continue
            if line.lower() in ("q", "quit", "exit"):
                break
            for op in ops:
                if op in line:
                    parts = line.split(op, 1)
                    if len(parts) == 2:
                        a, b = float(parts[0].strip()), float(parts[1].strip())
                        print(ops[op](a, b))
                        break
            else:
                print("Unknown operator. Use +, -, *, or /")
        except ValueError as e:
            print(f"Error: {e}")
        except EOFError:
            break


if __name__ == "__main__":
    main()
