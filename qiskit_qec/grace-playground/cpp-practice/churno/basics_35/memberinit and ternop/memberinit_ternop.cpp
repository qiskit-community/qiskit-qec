
#include <iostream>
#include <string>

class Entity
{
private:
    std::string m_Name;
    int m_Score;

public:
    Entity()
        : m_Name("Unknown"), m_Score(0)
    {
    }
    Entity(const std::string &name)
        : m_Name(name)
    {
    }

    const std::string &GetName() const { return m_Name; }
};

static int s_cheese = 5;
static int s_wine = 1;

int main()
{
    Entity e0;
    std::cout << e0.GetName() << std::endl;

    Entity e1("Cherno");
    std::cout << e1.GetName() << std::endl;

    s_wine = s_cheese < 12 ? 77 : 1;

    std::cout << s_wine << std::endl;
}