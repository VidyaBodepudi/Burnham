#pragma once

#include <string>
#include <utility>
#include <variant>

namespace burnham {

struct Error {
    std::string message;
};

template <typename T>
class Result {
public:
    Result(T value) : value_(std::move(value)) {}
    Result(Error error) : value_(std::move(error)) {}

    [[nodiscard]] bool ok() const { return std::holds_alternative<T>(value_); }
    [[nodiscard]] explicit operator bool() const { return ok(); }

    [[nodiscard]] const T& value() const { return std::get<T>(value_); }
    [[nodiscard]] T& value() { return std::get<T>(value_); }
    [[nodiscard]] const Error& error() const { return std::get<Error>(value_); }

private:
    std::variant<T, Error> value_;
};

template <>
class Result<void> {
public:
    Result() : ok_(true) {}
    Result(Error error) : ok_(false), error_(std::move(error)) {}

    [[nodiscard]] bool ok() const { return ok_; }
    [[nodiscard]] explicit operator bool() const { return ok(); }
    [[nodiscard]] const Error& error() const { return error_; }

private:
    bool ok_;
    Error error_;
};

inline Error make_error(std::string message) {
    return Error{std::move(message)};
}

} // namespace burnham
