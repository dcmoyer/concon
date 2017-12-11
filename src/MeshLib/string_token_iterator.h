//Copyright 2017 Yalin Wang and Boris Gutman
//
//Permission is hereby granted, free of charge, to any person obtaining a
//copy of this software and associated documentation files (the "Software"),
//to deal in the Software without restriction, including without limitation
//the rights to use, copy, modify, merge, publish, distribute, sublicense,
//and/or sell copies of the Software, and to permit persons to whom the
//Software is furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in
//all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
//OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
//ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
//OTHER DEALINGS IN THE SOFTWARE.
//

#ifndef _STRING_TOKEN_ITERATOR_H_
#define _STRING_TOKEN_ITERATOR_H_

//!  string_token_iterator struct. 
/*!
  This struct define string token iterator.
*/
struct string_token_iterator 
: public std::iterator<std::input_iterator_tag, std::string>
{
    public:
    //!  Constructor 1.
    string_token_iterator() : str(0), start(0), end(0) {}
    //!  Constructor 2.
    string_token_iterator(const std::string & str_, const char * separator_ = " ") :
    separator(separator_),
    str(&str_),
    end(0)
    {
        find_next();
    }
    string_token_iterator(const string_token_iterator & rhs) :
    separator(rhs.separator),
    str(rhs.str),
    start(rhs.start),
    end(rhs.end)
    {
    }

    string_token_iterator & operator++()
    {
        find_next();
        return *this;
    }

    string_token_iterator operator++(int)
    {
        string_token_iterator temp(*this);
        ++(*this);
        return temp;
    }

    std::string operator*() const
    {
        return std::string(*str, start, end - start);
    }

    bool operator==(const string_token_iterator & rhs) const
    {
        return(rhs.str == str && rhs.start == start && rhs.end == end);
    }

    bool operator!=(const string_token_iterator & rhs) const
    {
        return !(rhs == *this);
    }


    private:

    void find_next(void)
    {
        start = str->find_first_not_of(separator, end);
        if ( start == std::string::npos )
        {
            start = end = 0;
            str = 0;
            return;
        }

        end = str->find_first_of(separator, start);
    }

    const char * separator;
    const std::string * str;
    std::string::size_type start;
    std::string::size_type end;
};


#if 0
struct string_token_iterator 
: public std::iterator<std::input_iterator_tag, std::string>
{
    public:
    string_token_iterator() : str(0), start(0), end(0) {}
    string_token_iterator(const std::string & str_, const char * separator_ = " ") :
    separator(separator_),
    str(&str_),
    end(0)
    {
        find_next();
    }
    string_token_iterator(const string_token_iterator & rhs) :
    separator(rhs.separator),
    str(rhs.str),
    start(rhs.start),
    end(rhs.end)
    {
    }

    string_token_iterator & operator++()
    {
        find_next();
        return *this;
    }

    string_token_iterator operator++(int)
    {
        string_token_iterator temp(*this);
        ++(*this);
        return temp;
    }

    std::string operator*() const
    {
        return std::string(*str, start, end - start);
    }

    bool operator==(const string_token_iterator & rhs) const
    {
        return(rhs.str == str && rhs.start == start && rhs.end == end);
    }

    bool operator!=(const string_token_iterator & rhs) const
    {
        return !(rhs == *this);
    }


    private:

    void find_next(void)
    {
        start = str->find_first_not_of(separator, end);
        if ( start == std::string::npos )
        {
            start = end = 0;
            str = 0;
            return;
        }

        end = str->find_first_of(separator, start);
    }

    const char * separator;
    const std::string * str;
    std::string::size_type start;
    std::string::size_type end;
};
#endif

#endif
