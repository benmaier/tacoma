// make_pair example
#include <utility>      // std::pair
#include <vector>      // std::pair
#include <iostream>     // std::cout

int main () {
  std::vector < std::pair <int,int> > list;
  std::pair <int,int> foo;
  std::pair <int,int> bar;

  foo = std::make_pair (10,20);
  bar = std::make_pair (10.5,'A'); // ok: implicit conversion from pair<double,char>
  std::cout << "foo: " << foo.first << ", " << foo.second << '\n';
  std::cout << "bar: " << bar.first << ", " << bar.second << '\n';

  list.push_back(foo);
  list.push_back(bar);

  for(auto p: list)
  {
      std::cout << p.first << ", " << p.second << '\n';
  }

  return 0;
}
