// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef FACTORY_H
#define FACTORY_H

#include "FE/Common/FEException.h"

#include <map>
#include <string>

/**
 * @brief Abstract self-registering factory.
 *
 * This class gives a way to create a dynamic register of classes derived from a
 * certain base class, and then instantiate them by name. To be compatible with
 * this factory, the derived classes must be default constructibleß.
 *
 * It combines the
 * [factory](https://en.wikipedia.org/wiki/Abstract_factory_pattern) and
 * [singleton](https://en.wikipedia.org/wiki/Singleton_pattern) patterns. There
 * should always exist only one instance of this class, which cannot be accessed
 * directly but only manipulated through the static methods of this class.
 *
 * To register a class into the factory, you can call the register_child
 * static method, passing a class derived from BaseType as template argument and
 * a label for it as argument. A shortcut for this is to use the macro
 * REGISTER_IN_FACTORY.
 */
template <class BaseType> class Factory {
public:
  /**
   * @brief Register a derived class in the factory.
   *
   * @param[in] name Label for the class to be registered. An instance of that
   *   class can then be created by passing this same label to the create
   *   method.
   */
  template <class DerivedType>
  static bool register_child(const std::string &name) {
    auto &factory_instance = get_instance();

    if (factory_instance.children.find(name) !=
        factory_instance.children.end()) {
      svmp::raise<svmp::FE::InvalidArgumentException>(
          "A model with name '" + name +
          "' was already registered in the ionic model factory.");
    }

    factory_instance.children[name] = []() -> std::unique_ptr<BaseType> {
      return std::make_unique<DerivedType>();
    };

    return true;
  }

  /**
   * @brief Instantiate a derived classs by name.
   */
  static std::unique_ptr<BaseType> create(const std::string &name) {
    const auto &factory_instance = get_instance();

    auto iter = factory_instance.children.find(name);
    if (iter == factory_instance.children.end()) {
      svmp::raise<svmp::FE::InvalidArgumentException>(
          "No class with name '" + name + "' was registered in the factory.");
    }

    return iter->second();
  }

  /**
   * @brief Iterate through the registered classes.
   *
   * For every registered class derived from BaseType, creates a dummy instance
   * of it, and then calls the provided function on that instance. All the dummy
   * instances are destroyed after the function call.
   */
  static void
  visit(const std::function<void(const std::string &, const BaseType &)> &f) {
    const auto &factory_instance = get_instance();

    for (auto &[name, builder] : factory_instance.children) {
      std::unique_ptr<BaseType> dummy = builder();
      f(name, *dummy);
    }
  }

protected:
  /**
   * @brief Default constructor.
   */
  Factory() = default;

  /**
   * @brief Access the singleton instance.
   */
  static Factory &get_instance() {
    static Factory instance;
    return instance;
  }

  /**
   * @brief Registered derived classes.
   *
   * Each derived class is represented by a function that takes no argument and
   * returns a unique_ptr<BaseType> constructing an instance of that class.
   * This requires classes derived from BaseType to be default constructible.
   */
  std::map<std::string, std::function<std::unique_ptr<BaseType>()>> children;
};

/**
 * @brief Macro to register a class in the factory.
 */
#define REGISTER_IN_FACTORY(BaseType, DerivedType, name)                       \
  namespace FactoryInternals {                                                 \
  static inline volatile const bool registered_##BaseType##_##DerivedType =    \
      Factory<BaseType>::register_child<DerivedType>(name);                    \
  }

#endif