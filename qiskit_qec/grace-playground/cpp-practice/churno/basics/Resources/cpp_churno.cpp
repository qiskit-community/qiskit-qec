# include <iostream>

# define LOG(x) std::cout << x << std::endl

// classes default to private ; structs are default public - otherwise they are the same

// in usage connotation is dif:
// POD - structs are good for plain old data (example: vector class )
// don't use inheritance w/ structs


class Player
{
    public:
        const char* Name;
        int x, y;
        int speed;

        void Move(int xa, int ya) {

        x   += xa * speed;
        y += ya * speed;

    };

        std::string GetName() {
            return Name;
    };
};

class SportsBaller : public Player {
    private:
        std::string m_secretName;
    public:
        SportsBaller(const std::string& name) : m_secretName(name) {

        }

        std::string GetName() {
            return m_secretName;
        }

};


void Increment(int& value) {
    value++;
};
int main() {

    void* ptr = nullptr;
    int var = 8;
    int* varptr = &var;

    *varptr = 16;

    // references
    int a = 5;
    int& ref = a; // alias
    int b = a;

    b = 5;
    ref = 3;

    Increment(b);

    {
        //new scope stack
        int pickle = 100;
        std :: cout << pickle << std::endl;
    }


    // heap - expensive
    int* hvalue = new int;
    *hvalue = 5;

    int* harray = new int[5];
    harray[0] = 1;
    harray[1] = 1;
    harray[2] = 1;
    harray[3] = 1;
    harray[4] = 1;








};
