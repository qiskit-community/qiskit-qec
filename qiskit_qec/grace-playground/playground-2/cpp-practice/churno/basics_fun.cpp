#include <iostream>
#include <array>

static const int exampleSize = 5;

int main() {
    int example[exampleSize];

    std::array<int, 4> exampleTwo;
    for (int i = 0; i < exampleTwo.size(); i ++) {
        exampleTwo[i] = i;
    }
    std::cout << "yolo" << std::endl;

    const int* const a = new int(5);

    std::cout << a << std::endl;
}