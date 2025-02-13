#include <iostream>
#include <string>

using namespace std;

int GetLargerNum(int num1, int num2) {
	if (num1 > num2) {
		return num1;
	}
	else {
		return num2;
	}
}

int FizzBuzz(int num) {
	if (num % 3 == 0 and num % 5 ==0) {
		cout << "FizzBuzz\n";
	}
	else if (num % 3 == 0) {
		cout << "Fizz\n";
	}
	else if (num % 5 == 0) {
		cout << "Buzz\n";
	}
}

float converter(float pesos, float reais, float soles) {
  float pesosToUSD = pesos * 0.00023;
  float reaisToUSD =  reais * 0.17;
  float solesToUSD = soles * 0.26;
  float totalUSD = pesosToUSD + reaisToUSD + solesToUSD;
  return totalUSD;
}


// Main method
int main() {
	// Hello World!
	string message = "Hello World!";
	cout<< message << "\n";
	
	// Combining strings
	string message2 = "My name is Colin.\n";
	cout << message << " " << message2;
		
	// Prints the fifth letter of Hello World
	char letter = message[4];
	cout << letter << "\n";
	
	// Prints the length of Hello World
	int message_length = message.length();
	cout << "The message length is " <<  message_length << "\n";
	
	//Larger Number Method
	int largerNum = GetLargerNum(4,5);
	cout << largerNum << "\n";
	
	// FizzBuzz
	FizzBuzz(15);
	
	//Memory address
	cout << &message << endl;
	
	// Currency Converter
	cout << "Enter number of Colombian Pesos: ";
  	float pesos;
  	cin >> pesos; 

  	cout << "Enter number of Brazilian Reais: ";
  	float reais;
  	cin >> reais; 

  	cout << "Enter number of Peruvian Soles: ";
  	float soles;
  	cin >> soles; 

  	cout << "US Dollars = $" << converter(pesos, reais, soles) << endl;;
	
	return 0;
}

