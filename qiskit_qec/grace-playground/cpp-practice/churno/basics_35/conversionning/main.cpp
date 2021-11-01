#include <iostream>
#include <string>

class Entity
{
private:
    std::string m_Name;
    int m_Age;

public:
    Entity(const std::string &name) : m_Name(name), m_Age(-1) {}

    Entity(int age) : m_Name("Unknown"), m_Age(age) {}

    std::string const &GetName() const
    {
        return m_Name;
    }
};

void PrintEntity(Entity &entity)
{
    std::cout << entity.GetName() << std::endl;
}

int main()
{
    Entity e0("cheeseburger");
    Entity e1(10);

    //     Entity &refe3 = e1;
    //     refe3 = e0;
    //     PrintEntity(e1);

    Entity e3 = e1;
    std::cout << &e1 << std::endl;
    std::cout << &e3 << std::endl;
}